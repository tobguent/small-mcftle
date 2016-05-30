#pragma once
#ifndef _FLOW_INCLUDE_ONCE
#define _FLOW_INCLUDE_ONCE

#include "math.hpp"

// Contains the double gyre flow (standard parameters from Shadden:  http://mmae.iit.edu/shadden/LCS-tutorial/examples.html )
class Flow
{
public:
	Flow() : _BoundsMin(0,0,0), _BoundsMax(2,1,1) {}

	// Checks whether a coordinate is contained in the (spatial) flow domain
	bool Contains(const Vec3d& pos) const {
		return _BoundsMin.x() < pos.x() && pos.x() < _BoundsMax.x()
			&& _BoundsMin.y() < pos.y() && pos.y() < _BoundsMax.y()
			&& _BoundsMin.z() < pos.z() && pos.z() < _BoundsMax.z();
	}

	// Intersects a ray with the domain and returns entry and exit distance
	bool IntersectRay(const Ray& ray, Vec2d& t) const
	{
		t = Vec2d(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
		double tmin = (_BoundsMin.x() - ray.Origin.x()) / ray.Direction.x();
		double tmax = (_BoundsMax.x() - ray.Origin.x()) / ray.Direction.x();
		if (tmin > tmax) std::swap(tmin, tmax);
		double tymin = (_BoundsMin.y() - ray.Origin.y()) / ray.Direction.y();
		double tymax = (_BoundsMax.y() - ray.Origin.y()) / ray.Direction.y();
		if (tymin > tymax) std::swap(tymin, tymax);
		if ((tmin > tymax) || (tymin > tmax))
			return false;
		if (tymin > tmin)
			tmin = tymin;
		if (tymax < tmax)
			tmax = tymax;
		double tzmin = (_BoundsMin.z() - ray.Origin.z()) / ray.Direction.z();
		double tzmax = (_BoundsMax.z() - ray.Origin.z()) / ray.Direction.z();
		if (tzmin > tzmax) std::swap(tzmin, tzmax);
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

	// Samples the flow
	Vec3d Sample(const Vec3d& pos, const double& time, bool& indomain) const
	{
		indomain = Contains(pos);
		if (!indomain) return Vec3d(0,0,0);
		
		const double A = 0.1;
		return Vec3d(-M_PI * A * sin(M_PI * f(pos.x(),time)) * cos(M_PI * pos.y()),
					  M_PI * A * cos(M_PI * f(pos.x(),time)) * sin(M_PI * pos.y()) * dfdx(pos.x(),time),
					  0.0);
	}

	// Fourth-order Runge Kutta
	Vec3d StepRK4(const Vec3d& pos, double& time, double h, bool& indomain) const
	{
		Vec3d sample0 = Sample(pos, time, indomain);							if (!indomain) return pos;
		Vec3d sample1 = Sample(pos + sample0*0.5*h, time + 0.5*h, indomain);	if (!indomain) return pos;
		Vec3d sample2 = Sample(pos + sample1*0.5*h, time + 0.5*h, indomain);	if (!indomain) return pos;
		Vec3d sample3 = Sample(pos + sample2*h, time + h, indomain);			if (!indomain) return pos;
		time += h;
		return pos + (1.0/6.0) * (sample0 + 2.0*(sample1+sample2) + sample3) * h;   
	}

	// Computes the flow map
	Vec3d FlowMap(const Vec3d& x, double startTime, double endTime, double h, bool& indomain) const
	{
		// initialize
		Vec3d pos = x;
		double time = startTime;

		// sanity check
		indomain = Contains(x);								if (!indomain) return pos;

		// trace the particle
		while (time + h < endTime) {
			pos = StepRK4(pos, time, h, indomain);			if (!indomain) return pos;
		}
		pos = StepRK4(pos, time, endTime-time, indomain);	if (!indomain) return pos;
		return pos;
	}

	// Samples the FTLE value at a certain location x.
	double SampleFTLE(const Vec3d& x, double startTime, double endTime, double stepSize, const Vec3d& sepDistance) const
	{
		Vec3d pos = x;
		double t0 = startTime;
		double t1 = endTime;
		bool indomain = true;

		// compute finite differences
		Vec3d posX1 = FlowMap(pos + Vec3d(sepDistance.x(),0,0),t0,t1,stepSize,indomain);	if (!indomain) return 0;
		Vec3d posX0 = FlowMap(pos - Vec3d(sepDistance.x(),0,0),t0,t1,stepSize,indomain);	if (!indomain) return 0;

		Vec3d posY1 = FlowMap(pos + Vec3d(0,sepDistance.y(),0),t0,t1,stepSize,indomain);	if (!indomain) return 0;
		Vec3d posY0 = FlowMap(pos - Vec3d(0,sepDistance.y(),0),t0,t1,stepSize,indomain);	if (!indomain) return 0;

		Vec3d posZ1 = FlowMap(pos + Vec3d(0,0,sepDistance.z()),t0,t1,stepSize,indomain);	if (!indomain) return 0;
		Vec3d posZ0 = FlowMap(pos - Vec3d(0,0,sepDistance.z()),t0,t1,stepSize,indomain);	if (!indomain) return 0;

		Vec3d dx = (posX1 - posX0) / (2*sepDistance.x());
		Vec3d dy = (posY1 - posY0) / (2*sepDistance.y());
		Vec3d dz = (posZ1 - posZ0) / (2*sepDistance.z());
	
		// setup Jacobian and Cauchy Green tensor P
		Mat3x3d J = Mat3x3d(dx.x(), dy.x(), dz.x(),   dx.y(), dy.y(), dz.y(),  dx.z(), dy.z(), dz.z());
		Mat3x3d JT = J;	JT.transpose();
		Mat3x3d P = JT*J;

		// compute largest eigenvalue and finally the FTLE value
		double Lmax = P.LambdaMax();
		double ftle = 1 / std::abs(t1-t0) * log(sqrt(Lmax));
		return ftle;
	}

private:

	Vec3d _BoundsMin;
	Vec3d _BoundsMax;

	static double a(double t) {
		const double EPS = 0.25;
		const double OMEGA = 2*M_PI/10;
		return EPS * sin(OMEGA * t);
	}

	static double b(double t) {
		const double EPS = 0.25;
		const double OMEGA = 2*M_PI/10;
		return 1 - 2 * EPS * sin(OMEGA * t);
	}

	static double f(double x, double t) {
		return a(t)*x*x + b(t) * x;
	}

	static double dfdx(double x, double t) {
		return 2*a(t)*x + b(t);
	}
};


#endif