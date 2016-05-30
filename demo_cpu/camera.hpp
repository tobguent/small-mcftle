#pragma once
#ifndef _CAMERA_INCLUDE_ONCE
#define _CAMERA_INCLUDE_ONCE

#include "math.hpp"
#include <random>

// Pinhole camera.
class Camera
{
	public:
		// Creates a pinhole camera. Up vector should be normalized.
		Camera(const Vec3d& eye, const Vec3d& lookAt, const Vec3d& up, double fov, const Vec2i& screenResolution) : _Eye(eye), _ScreenResolution(screenResolution)
		{
			_W = lookAt - _Eye;
			double wlen = _W.length();
			_U = cross(_W, up).normalize();
			_V = cross(_U, _W).normalize();
			double ulen = wlen * tan(fov/180.0*M_PI);
			_U = _U * ulen;
			double aspect = (double)screenResolution.x/(double)screenResolution.y;
			double vlen = wlen * tan(fov/aspect/180.0*M_PI);
			_V = _V * vlen;
		}

		// Generates a ray through a random position in a pixel.
		Ray GenerateRay(const Vec2i& pxCoord) 
		{
			double jitterx = _RndDistribution(_RndGenerator);
			double jittery = _RndDistribution(_RndGenerator);
			double dx = ((pxCoord.x + 1) + jitterx) / ((double)_ScreenResolution.x + 1.0) * 2.0 - 1.0;
			double dy = ((pxCoord.y + 1) + jittery) / ((double)_ScreenResolution.y + 1.0) * 2.0 - 1.0;
			Vec3d ray_origin = _Eye;
			Vec3d ray_direction = (dx*_U + dy*_V + _W).normalize();
			return Ray(ray_origin, ray_direction);
		}

	private:
		std::default_random_engine _RndGenerator;
		std::uniform_real_distribution<double> _RndDistribution;

		Vec3d _Eye;					// position of the camera
		Vec3d _U;					// right vector
		Vec3d _V;					// up vector
		Vec3d _W;					// forward vector
		Vec2i _ScreenResolution;	// resolution of the viewport
};

#endif