#pragma once
#ifndef _MATH_INCLUDE_ONCE
#define _MATH_INCLUDE_ONCE

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <algorithm>

// Class that represents a 3D vector
class Vec3d
{
	public:
	Vec3d(const Vec3d& other) { *this = other; }
	Vec3d() { x() = y() = z() = 0; }
	explicit Vec3d(double scalar) { x() = y() = z() = scalar; }
	explicit Vec3d(double x, double y, double z)
	{
		mData[0] = x;
		mData[1] = y;
		mData[2] = z;
	}

	const double& x() const { return mData[0]; }
	const double& y() const { return mData[1]; }
	const double& z() const { return mData[2]; }

	double& x() { return mData[0]; }
	double& y() { return mData[1]; }
	double& z() { return mData[2]; }

	Vec3d operator+(const Vec3d& other) const {
		return Vec3d(x()+other.x(), y()+other.y(), z()+other.z());
	}
	Vec3d operator-(const Vec3d& other) const {
		return Vec3d(x()-other.x(), y()-other.y(), z()-other.z());
	}
	Vec3d operator+(double val) const {
		return Vec3d(x()+val, y()+val, z()+val);
	}
	Vec3d operator-(double val) const {
		return Vec3d(x()-val, y()-val, z()-val);
	}
	Vec3d operator*(double val) const {
		return Vec3d(x()*val, y()*val, z()*val);
	}
	Vec3d operator/(double val) const {
		return Vec3d(x()/val, y()/val, z()/val);
	}
	Vec3d& operator+=(const Vec3d& other) {
		*this = *this + other;
		return *this;
	}
	Vec3d& operator-=(const Vec3d& other) {
		*this = *this - other;
		return *this;
	}
	Vec3d& operator*=(double val) {
		*this = *this * val;
		return *this;
	}
	Vec3d operator-() {
		return Vec3d(-x(), -y(), -z());
	}
	
	double& operator[](unsigned i) { return mData[i]; }
	const double& operator[](unsigned i) const { return mData[i]; }
	double length() const { return sqrt(x()*x()+y()*y()+z()*z()); }
	double lengthSquared() const { return x()*x()+y()*y()+z()*z(); }
	const Vec3d& normalize() {
		double l = length();
		if (l) *this *= (double)(1.0/l);
		return *this;
	}

	protected:
	double mData[3];
};

inline const Vec3d operator*(double val, const Vec3d& v) {
	return v * val;
}
inline double dot(const Vec3d& v1, const Vec3d& v2) { return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z(); }
inline Vec3d cross(const Vec3d& v1, const Vec3d& v2) {
	Vec3d t;
	t.x() = v1.y()*v2.z() - v1.z()*v2.y();
	t.y() = v1.z()*v2.x() - v1.x()*v2.z();
	t.z() = v1.x()*v2.y() - v1.y()*v2.x();
	return t;
}
inline Vec3d lerp(const Vec3d& x, const Vec3d& y, double t) { return x+(y-x)*t; }
inline Vec3d max(double v1, const Vec3d& v2) { return Vec3d(std::max(v1,v2.x()), std::max(v1,v2.y()), std::max(v1,v2.z())); }
inline Vec3d pow(const Vec3d& v, double exp) { return Vec3d(pow(v.x(), exp), pow(v.y(), exp), pow(v.z(), exp)); }

inline double saturate(double t) { return std::min(std::max(0.0, t), 1.0); }
inline double sign(double x) { if (x==0) return 0; return x < 0 ? -1 : 1; }

double cbrt(double x) { return (x < 0.0 ? -pow(-x, 1.0/3.0) : 0.0); }

// normal form: x^3 + Ax^2 + Bx + C = 0
// adapted from:  http://read.pudn.com/downloads21/sourcecode/graph/71499/gems/Roots3And4.c__.htm   (Jochen Schwarze, in Graphics Gems 1990)
double maxRoot(double A, double B, double C)
{
	// substitute x = y - A/3 to eliminate quadric term: x^3 +px + q = 0  
	double sq_A = A * A;
	double p = 1.0/3 * (- 1.0/3 * sq_A + B);  
	double q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);  

	// use Cardano's formula
	double cb_p = p * p * p;  
	double D = q * q + cb_p;  

	if (D==0)
	{
		if (q==0)  // one triple solution
		{
			return - 1.0/3.0 * A;
		}
		else  // one single and one double solution
		{
			double u = cbrt(-q);
			return std::max(2*u, -u) - 1.0/3.0 * A;
		}
	}
	else
	{
		if (D < 0)  // Casus irreducibilis: three real solutions
		{
			double phi = 1.0/3 * acos(-q / sqrt(-cb_p));  
			double t = 2 * sqrt(-p);  
			return std::max(std::max(
				   t * cos(phi),
				 - t * cos(phi + 3.14159265359 / 3)),
				 - t * cos(phi - 3.14159265359 / 3)) - 1.0/3.0 * A;
		}
		else // one real solution
		{
			double sqrt_D = sqrt(D);  
			double u =  cbrt(sqrt_D - q);  
			double v = -cbrt(sqrt_D + q); 
			return u + v - 1.0/3.0 * A;
		}
	}
}

// Class that represents a 3D matrix
class Mat3x3d
{
public:
	Mat3x3d() { mVec[0] = Vec3d(1,0,0); mVec[1] = Vec3d(0,1,0); mVec[2] = Vec3d(0,0,1); }
	Mat3x3d(double e00, double e01, double e02, double e10, double e11, double e12, double e20, double e21, double e22)
	{
	  e(0,0) = e00; e(0,1) = e01; e(0,2) = e02; 
	  e(1,0) = e10; e(1,1) = e11; e(1,2) = e12; 
	  e(2,0) = e20; e(2,1) = e21; e(2,2) = e22;
	}

	const double& e(int i, int j) const { return mVec[j][i]; }
	double& e(int i, int j) { return mVec[j][i]; }

	Mat3x3d& transpose()
	{
		double tmp;
		for(int i=0; i<3; ++i)
			for(int j=i; j<3; ++j)
			{
				tmp = e(j,i);
				e(j,i) = e(i,j);
				e(i,j) = tmp;
			}
		return *this;
	}

	static Mat3x3d multiply(const Mat3x3d& p, const Mat3x3d& q)
	{
	  Mat3x3d out;
	  out.e(0,0) = q.e(0,0)*p.e(0,0) + q.e(1,0)*p.e(0,1) + q.e(2,0)*p.e(0,2);
	  out.e(0,1) = q.e(0,1)*p.e(0,0) + q.e(1,1)*p.e(0,1) + q.e(2,1)*p.e(0,2);
	  out.e(0,2) = q.e(0,2)*p.e(0,0) + q.e(1,2)*p.e(0,1) + q.e(2,2)*p.e(0,2);

	  out.e(1,0) = q.e(0,0)*p.e(1,0) + q.e(1,0)*p.e(1,1) + q.e(2,0)*p.e(1,2);
	  out.e(1,1) = q.e(0,1)*p.e(1,0) + q.e(1,1)*p.e(1,1) + q.e(2,1)*p.e(1,2);
	  out.e(1,2) = q.e(0,2)*p.e(1,0) + q.e(1,2)*p.e(1,1) + q.e(2,2)*p.e(1,2);

	  out.e(2,0) = q.e(0,0)*p.e(2,0) + q.e(1,0)*p.e(2,1) + q.e(2,0)*p.e(2,2);
	  out.e(2,1) = q.e(0,1)*p.e(2,0) + q.e(1,1)*p.e(2,1) + q.e(2,1)*p.e(2,2);
	  out.e(2,2) = q.e(0,2)*p.e(2,0) + q.e(1,2)*p.e(2,1) + q.e(2,2)*p.e(2,2);
	  return out;
	}

	inline Mat3x3d operator*(const Mat3x3d& q) { return multiply(*this, q); }

	// Computes the largest eigenvalue of a 3x3 matrix.
	double LambdaMax() const
	{
		double a = this->e(0,0);  double b = this->e(0,1);  double c = this->e(0,2);
		double d = this->e(1,0);  double e = this->e(1,1);  double f = this->e(1,2);
		double g = this->e(2,0);  double h = this->e(2,1);  double i = this->e(2,2);

		// determinant has the following polynomial x*x*x + a2*x*x + a1*x + a0 = 0 = det(P)
		double a2 = -(i+e+a);
		double a1 = -(-e*i-a*i+f*h+c*g-a*e+b*d);
		double a0 = -(a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g));

		return maxRoot(a2, a1, a0);
	}

private:
	Vec3d mVec[3];
};

struct Ray {
	Ray(const Vec3d& origin, const Vec3d& direction) : Origin(origin), Direction(direction) {}
	Vec3d Origin;
	Vec3d Direction;
};

struct Vec2d {
	Vec2d() : x(0), y(0) {}
	Vec2d(double X, double Y) : x(X), y(Y) {}
	double x;
	double y;
};

struct Vec2i {
	Vec2i(int X, int Y) : x(X), y(Y) {}
	int x;
	int y;
};

#endif