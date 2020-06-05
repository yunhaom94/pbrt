#pragma once
#include "core/ray.h"

// RayDifferential is a subclass of Ray that contains additional information about two
// auxiliary rays.These extra rays represent camera rays offset by one sample in the x and y
// direction from the main ray on the film plane.
class RayDifferential : public Ray
{
public:
	bool hasDifferentials;
	Eigen::Vector3d rxOrigin, ryOrigin, rxDirection, ryDirection;

public:
	RayDifferential() { hasDifferentials = false; }
	RayDifferential(const Ray& ray) : Ray(ray) {
		hasDifferentials = false;
	}
	void ScaleDifferentials(double s);

private:

};

