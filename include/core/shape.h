#pragma once

#include <Eigen/Core>

#include "core/pbrt.h"
#include "core/bounding_boxes.h"
#include "transform.h"



// this is the base shape interface
// very similiar to HW3 and HW4 in the idea
// of using polymorphism
class Shape
{

public:
	Transform s;
	const Transform *ObjectToWorld, *WorldToObject;
	bool reverseOrientation;
	bool transformSwapsHandedness;
	bool orientationIsAuthoritative;
	

public:
	Shape() {}

	Shape(const Transform* ObjectToWorld,
		const Transform* WorldToObject,
		bool reverseOrientation);

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

	// simple check for intersection, no additional infp
	virtual bool IntersectP(const Ray& ray,
		bool testAlphaTexture = true) const;

	// surface area
	virtual double Area() const = 0;



	~Shape() {}

private:

};

