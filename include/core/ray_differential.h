#pragma once
#include "core/ray.h"

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

