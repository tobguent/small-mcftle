#ifndef _HLSL_PARAMETERS_
#define _HLSL_PARAMETERS_

#define NUM_THREADS 256

#ifndef _WIN32
cbuffer cbParam : register( b0 )
{
	int2 screenResolution;				// viewport resolution
	int maxIterationsPerPixel;			// number of progressive iterations per pixel
	float startTime;					// FTLE parameters : start time of integration
	// --
	float endTime;						// FTLE parameters : end time of integration
	float stepSize;						// FTLE parameters : integration step size
	float separationDistance;			// FTLE parameters : finite difference epsilon
	float lightRadiance;				// radiance emitted from the light source  (L_e)
	// --
	float3 lightDirection;				// direction in which light emits photons (omega_L)
	float majorantExtinction;			// biggest possible extinction
	// --
	float2 boundsExtinction;			// transfer function bounds for the extinction coefficient sigma_t
	float2 boundsAlbedo;				// transfer function bounds for the scattering albedo c
	// --
	float3 U;							// camera parameters: right vector
	float contrast;						// tonemapping parameters 
	// --
	float3 V;							// camera parameters: up vector
	float brightness;					// tonemapping parameters 
	// --
	float3 W;							// camera parameters: forward vector
	float gamma;						// tonemapping parameters
	// --
	float3 eye;							// camera parameters: position of the camera
	float visibility_minimum;			// value of a visibilty failure (can be adjusted to get some ambient color in there for vis purposes)
}

cbuffer cbVarParam : register( b1 )
{
	int offset;							// thread offset in the pixel statistics buffer
}
#endif

#endif