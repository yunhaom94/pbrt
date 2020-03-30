#pragma once

#include "core/pbrt.h"
#include "core/shape.h"

class Sphere : public Shape
{
public:
	   const double radius;
	// part of the sphere can be cut off
	   const double zMin, zMax, thetaMin, thetaMax, phiMax;

public:

	Sphere(const Transform* ObjectToWorld,
		const Transform* WorldToObject,
		bool reverseOrientation,
		double radius,
		double zMin,
		double zMax,
		double phiMax);

	// return an bounding box around the shpere
	Bounds3d ObjectBound() const;

	// hw2
	bool Intersect(const Ray& r, double& tHit,
		SurfaceInteraction* isect, bool testAlphaTexture) const;

	// bool Sphere::IntersectP(const Ray &r, bool testAlphaTexture) const

	double Area() const;

	~Sphere() {}

private:

};

