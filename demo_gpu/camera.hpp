#ifndef _CAMERA_INCLUDE_ONCE
#define _CAMERA_INCLUDE_ONCE

#include "math.hpp"

// Initializes a pinhole camera.
void InitializeCamera(const Vec3f& eye, const Vec3f& lookAt, const Vec3f& up, float fov, const Vec2i& screenResolution, Vec3f& U, Vec3f& V, Vec3f& W, Vec3f& oEye)
{
	W = lookAt - eye;
	float wlen = length(W);
	U = normalize(cross(W, up));
	V = normalize(cross(U, W));
	float ulen = wlen * tanf(fov / 180.0f*(float)M_PI);
	U = U * ulen;
	float aspect = (float)screenResolution.x / (float)screenResolution.y;
	float vlen = wlen * tanf(fov / aspect / 180.0f*(float)M_PI);
	V = V * vlen;
	oEye = eye;
}

#endif