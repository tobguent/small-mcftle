#include "params.hlsli"
#include "flow.hlsli"
#include "random.hlsli"

// define resources
RWStructuredBuffer<float3> AccumRadiance : register( u0 );
RWStructuredBuffer<float> AccumCount : register( u1 );
RWStructuredBuffer<int2> RandomSeeds : register( u2 );

#define PI 3.14159265359

// ----------------------------------------------------------------------------------------------------------
// Transfer functions
// ----------------------------------------------------------------------------------------------------------
float ExtinctionCoefficient(float t) {
	return saturate((t-boundsExtinction.x)/(boundsExtinction.y-boundsExtinction.x)) * majorantExtinction;
}

float3 ScatteringAlbedo(float t) {
	t = 1-saturate((t-boundsAlbedo.x)/(boundsAlbedo.y-boundsAlbedo.x));
	float3 color = 
		(((float3(5.0048, 7.4158, 6.1246)*t + 
		float3(-8.0915, -15.9415, -16.2287))*t + 
		float3(1.1657, 7.4696, 11.9910))*t + 
		float3(1.4380, 1.2767, -1.4886))*t + 
		float3(0.6639, -0.0013, 0.1685);
	return color;
}

// ----------------------------------------------------------------------------------------------------------
// View ray generation and domain bounding box intersection
// ----------------------------------------------------------------------------------------------------------
void GenerateRay(int2 pxCoord, inout int2 seed, out float3 ray_origin, out float3 ray_direction)  {
	float jitterx = rnd(seed);
	float jittery = rnd(seed);
	float dx = ((pxCoord.x + 1) + jitterx) / ((float)screenResolution.x + 1.0) * 2.0 - 1.0;
	float dy = ((pxCoord.y + 1) + jittery) / ((float)screenResolution.y + 1.0) * 2.0 - 1.0;
	ray_origin = eye;
	ray_direction = normalize(dx*U + dy*V + W);
}

void swap(inout float a, inout float b) { float t = a; a = b; b = t; }

bool IntersectRay(float3 ray_origin, float3 ray_direction, inout float2 t)
{
	float3 _BoundsMin = float3(0,0,0);
	float3 _BoundsMax = float3(2,1,1);
	t = float2(-1000000, 1000000);
	float tmin = (_BoundsMin.x - ray_origin.x) / ray_direction.x;
	float tmax = (_BoundsMax.x - ray_origin.x) / ray_direction.x;
	if (tmin > tmax) swap(tmin, tmax);
	float tymin = (_BoundsMin.y - ray_origin.y) / ray_direction.y;
	float tymax = (_BoundsMax.y - ray_origin.y) / ray_direction.y;
	if (tymin > tymax) swap(tymin, tymax);
	if ((tmin > tymax) || (tymin > tmax))
		return false;
	if (tymin > tmin)
		tmin = tymin;
	if (tymax < tmax)
		tmax = tymax;
	float tzmin = (_BoundsMin.z - ray_origin.z) / ray_direction.z;
	float tzmax = (_BoundsMax.z - ray_origin.z) / ray_direction.z;
	if (tzmin > tzmax) swap(tzmin, tzmax);
	if ((tmin > tzmax) || (tzmin > tmax))
		return false;
	if (tzmin > tmin)
		tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;
	if ((tmin > t.y) || (tmax < t.x)) return false;
	if (t.x < tmin) t.x = tmin;
	if (t.y > tmax) t.y = tmax;
	return true;
}

// ----------------------------------------------------------------------------------------------------------
// Free Path Sampling
// ----------------------------------------------------------------------------------------------------------
struct FreePathResult {
	float d;	// free flight distance until the scattering event
	float dmax;	// maximal flight distance (if d<=dmax, we have a scattering event)
	float ftle;	// scalar field value at scattering event
};

// performs a free path sampling, aka. Woodcock tracking
inline FreePathResult FreePathSampling(float3 ray_origin, float3 ray_direction, inout int2 seed)
{
	// initialize result data structure to 'no hit' (d > dmax).
	FreePathResult res;
	res.d = 3.402823466e+38F;
	res.dmax = 0;
	res.ftle = 0;

	// first, test if the ray hits the domain at all (hits are sorted, thus if global_hit.y<0 then both are negative)
	float2 global_hit;
	if (!IntersectRay(ray_origin, ray_direction, global_hit) || global_hit.y < 0) 
		return res;

	// set extinction coefficient (global upper bound) and entry and exit of the run
	float km = majorantExtinction;
	float dmin = max(0.0,global_hit.x);
	res.dmax = global_hit.y;

	// perform the free path run (Alg. 1.)
	res.d = dmin - log(1-rnd(seed)) / km;
	[allow_uav_condition]
	while (res.d <= res.dmax && ExtinctionCoefficient(res.ftle = sampleFTLE(ray_origin + res.d * ray_direction)) / km < rnd(seed))
		res.d = res.d - log(1-rnd(seed)) / km;
	return res;
}

// ----------------------------------------------------------------------------------------------------------
// Compute Shader
// ----------------------------------------------------------------------------------------------------------
[numthreads(NUM_THREADS, 1, 1)]
void main( uint DTid : SV_DispatchThreadID, uint Tid : SV_GroupThreadID, uint Gid : SV_GroupID )
{
	uint px = DTid+offset;
	if (px >= (uint)(screenResolution.x*screenResolution.y)) return;

	int2 seed = RandomSeeds[px];
	
	// select pixel and generate a ray
	int2 pxCoord = int2(px%screenResolution.x, px/screenResolution.x);
	float3 viewRay_origin, viewRay_direction;
	GenerateRay(pxCoord, seed, viewRay_origin, viewRay_direction);
	
	// intersect ray with domain (entry and exit)
	float2 t = float2(0,0);
	if (!IntersectRay(viewRay_origin, viewRay_direction, t)) return;

	// free path sampling to find a scattering event
	FreePathResult fprView = FreePathSampling(viewRay_origin, viewRay_direction, seed);

	float3 radiance = float3(1,1,1);	// initialize with radiance of the background (white)

	// if we found an event
	if (fprView.d <= fprView.dmax)
	{
		// compute the scattering albedo
		float3 albedo = ScatteringAlbedo(fprView.ftle);

		// compute scattering location x_s and generate shadow ray
		float3 x_s = viewRay_origin + fprView.d * viewRay_direction;
		if (!IntersectRay(x_s, -lightDirection, t)) return;
			
		// free path sampling to find a scattering event on the shadow ray
		FreePathResult fprLight = FreePathSampling(x_s, -lightDirection, seed);

		// if we have scattered, we are in shadow. (actually visibility would be zero now, but for vis purposes we keep a little bit.)
		float visibility = 1;
		if (fprLight.d <= fprLight.dmax)
			visibility = visibility_minimum;

		// phase function
		const float phase = 1.0/(4.0*PI);		// isotropic phase function
			
		// compute in-scattered radiance
		radiance = albedo * phase * visibility * lightRadiance;
	}
	// add the result to the pixel statistics
	AccumRadiance[px] += radiance;
	AccumCount[px] += 1;
			
	// Store random seed
	RandomSeeds[px] = seed;
}

