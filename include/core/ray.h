#pragma once

#include "core/pbrt.h"

class Ray
{
public:
	// parametric ray r(t) = o + td
	// default to (0,0,0)
	Point3f o = Point3f(0, 0, 0);
	Vector3f d = Vector3f(0, 0, 0);
	// mutable means is able to change 
	// even if it marked as const
	mutable Float tMax;
	// for animation
	Float time;
	const Medium *medium;

public:
	Ray() : tMax(Infinity), time(0), medium(nullptr) {}

	Ray(const Point3f& o, 
		const Vector3f& d, 
		Float tMax = Infinity,
		Float time = 0, 
		const Medium* medium = nullptr)
		: o(o), d(d), tMax(tMax), time(time), medium(medium) {}

	~Ray() {}

	// overload () operator
	Point3f operator()(Float t) const;

};


// RayDifferential is a subclass of Ray that contains additional information about two
// auxiliary rays. These extra rays represent camera rays offset by one sample in the x and y
// direction from the main ray on the film plane.
// Used to texture anti-aliasing from averaging the area of texture covered by 3 rays.
class RayDifferential : public Ray
{
public:
	bool hasDifferentials;
	Point3f rxOrigin, ryOrigin;
	Vector3f rxDirection, ryDirection;

public:
	RayDifferential() { hasDifferentials = false; }
	RayDifferential(const Ray& ray) : Ray(ray) 
	{
		hasDifferentials = false;
	}
	RayDifferential(const Point3f& o, const Vector3f& d, Float tMax = Infinity,
		Float time = 0.f, const Medium* medium = nullptr)
		: Ray(o, d, tMax, time, medium) {
		hasDifferentials = false;
	}
	// update differential rays for an estimated sample spacing of s.
	void ScaleDifferentials(Float s);

};



