#pragma once

#include "core/pbrt.h"
#include "core/shape.h"

class Sphere : public Shape
{
public:
	const Float radius;
	// if over the section of sphere is cut off
	const Float zMin, zMax;
	const Float thetaMin, thetaMax, phiMax;

public:

	Sphere(const Transform* ObjectToWorld,
		const Transform* WorldToObject,
		bool reverseOrientation,
		Float radius,
		Float zMin,
		Float zMax,
		Float phiMax);

	~Sphere() {}

	// return an bounding box around the shpere
	Bounds3f ObjectBound() const;

	bool Intersect(const Ray& r, Float* tHit,
		SurfaceInteraction* isect, bool testAlphaTexture) const;

	// TODO: try this later, intersect with 418 hw2
	//bool Intersect418(const Ray& r, double& tHit, SurfaceInteraction* isect, bool testAlphaTexture) const; 


	bool IntersectP(const Ray& r, bool testAlphaTexture) const;

	double Area() const;

	Interaction Sphere::Sample(const Point2f& u) const;

	Interaction Sample(const Interaction& ref, const Point2f& u) const;

	Float Pdf(const Interaction& ref, const Vector3f& wi) const;

private:

};

