#pragma once

#include "core/pbrt.h"
#include "core/ray.h"
#include "core/ray_differential.h"

// holds all of the sample values needed to specify a camera ray.
/*
struct CameraSample
{
	Eigen::Vector2d pFilm;
	Eigen::Vector2d  pLens;
	double time;
};

class Camera
{

public:
	AnimatedTransform CameraToWorld;
	const double shutterOpen, shutterClose;
	Film* film;
	const Medium* medium;

public:
	Camera(const AnimatedTransform& CameraToWorld, double shutterOpen,
		double shutterClose, Film* film, const Medium* medium);

	virtual double GenerateRay(const CameraSample& sample,
		Ray* ray) const = 0;

	// computes a main ray like GenerateRay() but
	// also computes the corresponding rays for pixels shifted one pixel in the xand y directions
	// on the film plane
	double Camera::GenerateRayDifferential(const CameraSample& sample,
		RayDifferential* rd) const;

	~Camera();

private:

};

*/