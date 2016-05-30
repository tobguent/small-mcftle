#pragma once
#ifndef _MATH_INCLUDE_ONCE
#define _MATH_INCLUDE_ONCE

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <algorithm>

struct Vec3f {
	Vec3f() : x(0), y(0), z(0) {}
	Vec3f(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
	float x;
	float y;
	float z;
};

Vec3f operator+(const Vec3f& a, const Vec3f& b) { return Vec3f(a.x + b.x, a.y + b.y, a.z + b.z); }
Vec3f operator-(const Vec3f& a, const Vec3f& b) { return Vec3f(a.x - b.x, a.y - b.y, a.z - b.z); }

inline Vec3f operator+(const Vec3f& vec, float val) { return Vec3f(vec.x + val, vec.y + val, vec.z + val); }
inline Vec3f operator-(const Vec3f& vec, float val) { return Vec3f(vec.x - val, vec.y - val, vec.z - val); }
inline Vec3f operator*(const Vec3f& vec, float val) { return Vec3f(vec.x * val, vec.y * val, vec.z * val); }
inline Vec3f operator/(const Vec3f& vec, float val) { return Vec3f(vec.x / val, vec.y / val, vec.z / val); }
inline Vec3f max(float v1, const Vec3f& v2) { return Vec3f(std::max(v1, v2.x), std::max(v1, v2.y), std::max(v1, v2.z)); }
inline Vec3f pow(const Vec3f& v, float exp) { return Vec3f(pow(v.x, exp), pow(v.y, exp), pow(v.z, exp)); }

inline float length(const Vec3f& v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
inline Vec3f normalize(const Vec3f& v) {
	float l = length(v);
	return l!=0.f ? v / l : v;
}
inline Vec3f cross(const Vec3f& v1, const Vec3f& v2) {
	return Vec3f(v1.y*v2.z - v1.z*v2.y,
		v1.z*v2.x - v1.x*v2.z,
		v1.x*v2.y - v1.y*v2.x);
}

struct Vec2f {
	Vec2f() : x(0), y(0) {}
	Vec2f(float X, float Y) : x(X), y(Y) {}
	float x;
	float y;
};

struct Vec2i {
	Vec2i(int X, int Y) : x(X), y(Y) {}
	int x;
	int y;
};

#endif