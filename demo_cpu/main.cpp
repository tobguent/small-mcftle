#include <omp.h>
#include "camera.hpp"
#include "frame.hpp"
#include "flow.hpp"
#include "timer.hpp"

// =============================================================
// set parameters
// =============================================================
Vec2i screenResolution(200, 100);											// viewport resolution
int maxIterationsPerPixel(100);												// number of progressive iterations per pixel
double startTime(0), endTime(40), stepSize(0.1), separationDistance(1e-6);	// FTLE parameters
Vec3d lightDirection(-1/sqrt(3),-1/sqrt(3),-1/sqrt(3));						// normalized direction in which light emits photons (omega_L)
double lightRadiance(12);													// radiance emitted from the light source  (L_e)
Vec2d boundsExtinction(0.05, 0.235);										// transfer function bounds for the extinction coefficient sigma_t
Vec2d boundsAlbedo(0.05, 0.235);											// transfer function bounds for the scattering albedo c
double majorantExtinction(100);												// biggest possible extinction
double contrast(1.85), brightness(0.25), gamma_val(1.4);					// tonemapping parameters
double visibility_minimum = 0.0;											// value of a visibilty failure (can be adjusted to get some ambient color in there for vis purposes)

std::default_random_engine rng;
std::uniform_real_distribution<double> rnd(0, 1);

// =============================================================
// transfer functions
// =============================================================
double ExtinctionCoefficient(double t) {
	return saturate((t-boundsExtinction.x)/(boundsExtinction.y-boundsExtinction.x)) * majorantExtinction;
}
Vec3d ScatteringAlbedo(double t) {
	t = 1-saturate((t-boundsAlbedo.x)/(boundsAlbedo.y-boundsAlbedo.x));	
	Vec3d color = 
		(((Vec3d(5.0048, 7.4158, 6.1246)*t +
		Vec3d(-8.0915, -15.9415, -16.2287))*t + 
		Vec3d(1.1657, 7.4696, 11.9910))*t + 
		Vec3d(1.4380, 1.2767, -1.4886))*t + 
		Vec3d(0.6639, -0.0013, 0.1685);
	return color;
}

// =============================================================
// free path sampling
// =============================================================
struct FreePathResult {
	double d;		// free flight distance until the scattering event occured
	double dmax;	// maximal flight distance (if d<=dmax, we have a scattering event)
	double ftle;	// scalar field value at scattering event
};

// performs a free path sampling, aka. Woodcock tracking
FreePathResult FreePathSampling(const Ray& ray, Flow& flow)
{
	// initialize result data structure to 'no hit' (d > dmax).
	FreePathResult res;
	res.d = std::numeric_limits<double>::max();
	res.dmax = 0;
	res.ftle = 0;

	// first, test if the ray hits the domain at all (hits are sorted, thus if global_hit.y<0 then both are negative)
	Vec2d global_hit;
	if (!flow.IntersectRay(ray, global_hit) || global_hit.y < 0) 
		return res;

	// set extinction coefficient (global upper bound) and entry and exit of the run
	double km = majorantExtinction;
	double dmin = std::max(0.0,global_hit.x);
	res.dmax = global_hit.y;

	// perform the free path run (Alg. 1)
	res.d = dmin - log(1-rnd(rng)) / km;
	while (res.d <= res.dmax && ExtinctionCoefficient(res.ftle = flow.SampleFTLE(ray.Origin + res.d * ray.Direction, startTime, endTime, stepSize, Vec3d(separationDistance))) / km < rnd(rng))
		res.d = res.d - log(1-rnd(rng)) / km;
	return res;
}

// =============================================================
// =============================================================
// =============================================================
int main()
{
	printf("Small MCFTLE: Computing...\n");

	// =============================================================
	// initialize
	// =============================================================
	Vec3d eye(1,0.5,2.3), lookAt(1,0.5,0.5), up(0,1,0);		// camera parameters
	double fov = 45.0;
	Camera camera(eye, lookAt, up, fov, screenResolution);	// pinhole camera
	Frame frame(screenResolution);							// stores pixel statistics
	Flow flow;												// vector field
	Timer timer;											// timer to measure runtime

	// =============================================================
	// progressive rendering
	// =============================================================
	for (int it=0; it<maxIterationsPerPixel; ++it)
	{
		timer.tic();	// start runtime measurement

		#ifdef NDEBUG
		#pragma omp parallel for schedule(dynamic,16)
		#endif
		for (int px=0; px<screenResolution.x*screenResolution.y; ++px)
		{
			// select pixel and generate a ray
			Vec2i pxCoord(px%screenResolution.x, px/screenResolution.x);
			Ray viewRay = camera.GenerateRay(pxCoord);

			// intersect ray with domain (entry and exit)
			Vec2d t;
			if (!flow.IntersectRay(viewRay, t))
				continue;

			// free path sampling to find a scattering event
			FreePathResult fprView = FreePathSampling(viewRay, flow);
			
			Vec3d radiance(1,1,1);	// initialize with radiance of the background (white)

			// if we found a scattering event
			if (fprView.d <= fprView.dmax)
			{
				// compute the scattering albedo c
				Vec3d albedo = ScatteringAlbedo(fprView.ftle);

				// compute scattering location x_s and generate shadow ray
				Vec3d x_s = viewRay.Origin + fprView.d * viewRay.Direction;
				Ray shadowRay(x_s, -lightDirection);

				// intersect shadow ray with domain (determine exit)
				if (!flow.IntersectRay(shadowRay, t))
					continue;

				// free path sampling to find a scattering event on the shadow ray
				FreePathResult fprLight = FreePathSampling(shadowRay, flow);

				// if we have scattered, we are in shadow. (actually visibility would be zero now, but for vis purposes we keep a little bit.)
				double visibility = 1;
				if (fprLight.d <= fprLight.dmax)
					visibility = visibility_minimum;

				// phase function
				double phase = 1.0 / (4.0 * M_PI);		// isotropic phase function

				// compute in-scattered radiance
				radiance = albedo * phase * visibility * lightRadiance;
			}
			// add the result to the pixel statistics
			frame.AddPixel(pxCoord, radiance);
		}
		
		timer.toc();	// finish runtime measurement (prints result to console)

		// write intermediate result to file
		frame.ExportPfm("result.pfm");
		frame.ExportBmp("result.bmp", contrast, brightness, gamma_val);
		printf("Iteration: %i / %i\n", it+1, maxIterationsPerPixel);
	}
	printf("\nComplete!\n");
	return 0;
}