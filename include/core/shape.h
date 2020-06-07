#pragma once

#include "core/pbrt.h"


// This is the base shape interface used by Primitive to contain geometric info.
// All spheres are defined in a coordinate system where the center of
// the sphere is at the origin.
class Shape
{

public:
	const Transform *ObjectToWorld, *WorldToObject;
	bool reverseOrientation;
	bool transformSwapsHandedness;
	bool orientationIsAuthoritative;
	

public:
	Shape() {}

	Shape(const Transform* ObjectToWorld,
		const Transform* WorldToObject,
		bool reverseOrientation);

	~Shape() {}

	// return a bounding box around the object
	virtual Bounds3f ObjectBound() const = 0;

	// bounding box, but in world space
	Bounds3f WorldBound() const;

	// do ray-object intersection
	// in: ray
	// out: tHit
	//		isect - surface interaction point info
	//		testAlphaTexture - Some shape implementations support 
	//			cutting away some of their surfaces using a texture
	virtual bool Intersect(const Ray& ray, double& tHit,
		SurfaceInteraction* isect, bool testAlphaTexture = true) const = 0;

	// simple check for intersection, no additional info
	virtual bool IntersectP(const Ray& ray,
		bool testAlphaTexture = true) const;

	// surface area
	virtual Float Area() const = 0;


};

