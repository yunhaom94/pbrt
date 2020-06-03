#pragma once

#include "core/pbrt.h"


// TODO: 356
//records the position on the film for which the camera should
//generate the corresponding ray
struct CameraSample
{
	Point2d pFilm;
	Point2d pLens;
	double time;
};

class Camera
{

public:
	//AnimatedTransform CameraToWorld;
	//const double shutterOpen, shutterClose;
	Film* film;
	//const Medium* medium;


public:

	Camera() {}
	~Camera() {}


	Camera(const AnimatedTransform& CameraToWorld, double shutterOpen,
		double shutterClose, Film* film, const Medium* medium);

	virtual double GenerateRay(const CameraSample& sample,
		Ray* ray) const = 0;

	// computes a main ray like GenerateRay() but
	// also computes the corresponding rays for pixels shifted one pixel in the xand y directions
	// on the film plane
	double Camera::GenerateRayDifferential(const CameraSample& sample,
		RayDifferential* rd) const;


private:

};

