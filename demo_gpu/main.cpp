#include "math.hpp"
#include "d3d.hpp"
#include "frame.hpp"
#include "cbuffer.hpp"
#include "camera.hpp"
#include "params.hlsli"
#include "timer.hpp"

// =============================================================
// set parameters
// =============================================================

struct Params	// global constant parameters, sorted for 16 byte alignment
{
	Params() : 
		screenResolution(400, 200),
		maxIterationsPerPixel(10000),
		startTime(0), endTime(40), stepSize(0.1f), separationDistance(1e-6f),
		lightRadiance(12), lightDirection(-1/sqrtf(3), -1/sqrtf(3), -1/sqrtf(3)),
		boundsExtinction(0.05f, 0.235f),
		boundsAlbedo(0.05f, 0.235f),
		majorantExtinction(100),
		contrast(1.85f), brightness(0.25f), gamma(1.4f),
		visibility_minimum(0.0f)
	{}

	Vec2i screenResolution;				// viewport resolution
	int maxIterationsPerPixel;			// number of progressive iterations per pixel
	float startTime;					// FTLE parameters : start time of integration
	// --
	float endTime;						// FTLE parameters : end time of integration
	float stepSize;						// FTLE parameters : integration step size
	float separationDistance;			// FTLE parameters : finite difference epsilon
	float lightRadiance;				// radiance emitted from the light source  (L_e)
	// --
	Vec3f lightDirection;				// direction in which light emits photons (omega_L)
	float majorantExtinction;			// biggest possible extinction
	// --
	Vec2f boundsExtinction;				// transfer function bounds for the extinction coefficient sigma_t
	Vec2f boundsAlbedo;					// transfer function bounds for the scattering albedo c
	// --
	Vec3f U;							// camera parameters: right vector
	float contrast;						// tonemapping parameters 
	// --
	Vec3f V;							// camera parameters: up vector
	float brightness;					// tonemapping parameters 
	// --
	Vec3f W;							// camera parameters: forward vector
	float gamma;						// tonemapping parameters
	// --
	Vec3f eye;							// camera parameters: position of the camera
	float visibility_minimum;			// value of a visibilty failure (can be adjusted to get some ambient color in there for vis purposes)
};

struct VaryingParams  // params that change every dispatch call
{
	VaryingParams() : offset(0) {}
	int offset;							// thread offset in the pixel statistics buffer
};

// =============================================================
// =============================================================
// =============================================================
int main()
{
	printf("Small MCFTLE: Computing...\n");

	// =============================================================
	// initialize
	// =============================================================

	// Initialize Direct Compute
	D3D d3d;
	ID3D11Device* device = d3d.GetDevice();
	ID3D11DeviceContext* immediateContext = d3d.GetImmediateContext();

	// Initialize parameters
	ConstantBuffer<Params> params;
	Vec3f eye(1, 0.5f, 2.3f), lookAt(1, 0.5f, 0.5f), up(0, 1, 0);
	float fov = 45.0f;
	InitializeCamera(eye, lookAt, up, fov, params.Data.screenResolution, params.Data.U, params.Data.V, params.Data.W, params.Data.eye);
	params.Create(device);

	ConstantBuffer<VaryingParams> varparams;	// parameters that change every dispatch call
	varparams.Create(device);

	// Initialize frame for pixel statistics
	Frame frame(params.Data.screenResolution, device);

	// create timer to measure runtime
	Timer timer;

	// Load compute shader
	ID3D11ComputeShader* computeShader = NULL;
	d3d.LoadComputeShaderFromFile("mcftle.cso", &computeShader);

	// =============================================================
	// progressive rendering
	// =============================================================

	// Bind buffers
	ID3D11Buffer* cbs[] = { params.GetBuffer(), varparams.GetBuffer() };
	immediateContext->CSSetConstantBuffers(0, 2, cbs);
	immediateContext->CSSetShader(computeShader, NULL, 0);
	ID3D11UnorderedAccessView* uavs[] = { frame.GetUavAccumRadiance(), frame.GetUavAccumCount(), frame.GetUavRandomSeeds() };
	UINT initialCounts[] = { 0, 0, 0 };
	immediateContext->CSSetUnorderedAccessViews(0, 3, uavs, initialCounts);

	// Calculate dispatch size (breakup the computation of a frame into multiple dispatches)
	const int numThreadsPerGroup = NUM_THREADS;
	const int numThreadGroups = 1024;
	const int numPixels = params.Data.screenResolution.x * params.Data.screenResolution.y;
	int numSamplesPerDispatch = std::min(numThreadGroups * numThreadsPerGroup, numPixels);
	int iterationsPerLoop = numPixels / numSamplesPerDispatch;
	if (numPixels % numSamplesPerDispatch != 0) iterationsPerLoop++;
	int threadGroups = numSamplesPerDispatch / numThreadsPerGroup;
	if (numSamplesPerDispatch % numThreadsPerGroup != 0) threadGroups++;

	// Compute loop
	for (int it = 0; it < params.Data.maxIterationsPerPixel; ++it)
	{
		timer.tic();	// start runtime measurement

		varparams.Data.offset = 0;
		for (int el = 0; el < iterationsPerLoop; ++el)
		{
			varparams.UpdateBuffer(immediateContext);
			immediateContext->Dispatch(threadGroups, 1, 1);
			varparams.Data.offset += numSamplesPerDispatch;
		}
		
		// write intermediate result to file
		frame.CopyToHost(immediateContext);

		timer.toc();	// finish runtime measurement (prints result to console)

		frame.ExportPfm("result.pfm");
		frame.ExportBmp("result.bmp", params.Data.contrast, params.Data.brightness, params.Data.gamma);
		printf("Iteration: %i / %i\n", it + 1, params.Data.maxIterationsPerPixel);
	}

	if (computeShader) computeShader->Release();
	printf("\nComplete!\n");
	return 0;
}