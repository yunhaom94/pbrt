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
	virtual Bounds3d ObjectBound() const = 0;

	// bounding box, but in world space
	Bounds3d WorldBound() const;

	~Shape() {}

private:

};

