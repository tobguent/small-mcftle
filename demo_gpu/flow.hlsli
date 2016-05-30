#ifndef _HLSL_FTLE_3D_
#define _HLSL_FTLE_3D_

#include "params.hlsli"

// doube gyre flow
float a(float t) {
	const float PI = 3.14159265359;
	const float EPS = 0.25;
	const float OMEGA = 2*PI/10;
	return EPS * sin(OMEGA * t);
}

float b(float t) {
	const float PI = 3.14159265359;
	const float EPS = 0.25;
	const float OMEGA = 2*PI/10;
	return 1 - 2 * EPS * sin(OMEGA * t);
}

float f(float x, float t) {
	return a(t)*x*x + b(t) * x;
}

float dfdx(float x, float t) {
	return 2*a(t)*x + b(t);
}

float3 sampleFlow(float3 pos, float t)
{
	const float A = 0.1;
	const float PI = 3.14159265359;

	float3 v = float3(
		-PI * A * sin(PI * f(pos.x,t)) * cos(PI * pos.y),
		 PI * A * cos(PI * f(pos.x,t)) * sin(PI * pos.y) * dfdx(pos.x,t),
		 0.0);

	return v;
}

float3 stepRK4(float3 pos, inout float time, float h) 
{
	float3 sample0 = sampleFlow(pos,time);
	float3 sample1 = sampleFlow(pos+sample0 * 0.5*h, time+0.5*h);
	float3 sample2 = sampleFlow(pos+sample1 * 0.5*h, time+0.5*h);
	float3 sample3 = sampleFlow(pos+sample2 * h, time+h);
	time += h;
	return pos + (1.0/6.0) * (sample0 + 2*(sample1+sample2) + sample3) * h;
}

float3 flowmap(float3 pos, float time, float endTime)
{
	[allow_uav_condition]
	while (time + stepSize < endTime)
		pos = stepRK4(pos, time, stepSize);
	return stepRK4(pos, time, endTime-time);
}

float cbrt(float x) { return (x < 0.0 ? -pow(abs(-x), 1.0/3.0) : 0.0); }

// normal form: x^3 + Ax^2 + Bx + C = 0
// adapted from:  http://read.pudn.com/downloads21/sourcecode/graph/71499/gems/Roots3And4.c__.htm   (Jochen Schwarze, in Graphics Gems 1990)
float maxRoot(float A, float B, float C)
{
	// substitute x = y - A/3 to eliminate quadric term: x^3 +px + q = 0  
	float sq_A = A * A;
	float p = 1.0/3 * (- 1.0/3 * sq_A + B);  
	float q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);  

	// use Cardano's formula
	float cb_p = p * p * p;  
	float D = q * q + cb_p;  

	if (D==0)
	{
		if (q==0)  // one triple solution
		{
			return - 1.0/3.0 * A;
		}
		else  // one single and one double solution
		{
			float u = cbrt(-q);
			return max(2*u, -u) - 1.0/3.0 * A;
		}
	}
	else
	{
		if (D < 0)  // Casus irreducibilis: three real solutions
		{
			float phi = 1.0/3 * acos(-q / sqrt(-cb_p));  
			float t = 2 * sqrt(-p);  
			return max(max(
				   t * cos(phi),
				 - t * cos(phi + 3.14159265359 / 3)),
				 - t * cos(phi - 3.14159265359 / 3)) - 1.0/3.0 * A;
		}
		else // one real solution
		{
			float sqrt_D = sqrt(D);  
			float u =  cbrt(sqrt_D - q);  
			float v = -cbrt(sqrt_D + q); 
			return u + v - 1.0/3.0 * A;
		}
	}
}

float sampleFTLE(float3 pos)
{
	float t0 = startTime;
	float t1 = endTime;
	float hh = separationDistance;

	float3 dx = (flowmap(pos + float3(hh,0,0),t0,t1) - flowmap(pos - float3(hh,0,0),t0,t1)) / (2*hh);
	float3 dy = (flowmap(pos + float3(0,hh,0),t0,t1) - flowmap(pos - float3(0,hh,0),t0,t1)) / (2*hh);
	float3 dz = (flowmap(pos + float3(0,0,hh),t0,t1) - flowmap(pos - float3(0,0,hh),t0,t1)) / (2*hh);

	float3x3 J = float3x3(dx.x, dy.x, dz.x,   dx.y, dy.y, dz.y,  dx.z, dy.z, dz.z);
	float3x3 P = mul(transpose(J),J);

	float a = P._11;  float b = P._12;  float c = P._13;
	float d = P._21;  float e = P._22;  float f = P._23;
	float g = P._31;  float h = P._32;  float i = P._33;

	// determinant has the following polynomial x*x*x + a2*x*x + a1*x + a0 = 0 = det(P)
	float a2 = -(i+e+a);
	float a1 = -(-e*i-a*i+f*h+c*g-a*e+b*d);
	float a0 = -(a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g));

	float Lmax = maxRoot(a2, a1, a0);
	float ftle = 1 / abs(t1-t0) * log(sqrt(Lmax));
	return ftle;
}
#endif