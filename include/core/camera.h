#pragma once

#include "core/pbrt.h"
#include "core/transform.h"


//records the position on the film for which the camera should
//generate the corresponding ray
struct CameraSample
{
	//the point on the film to which the generated ray carries radiance
	Point2f pFilm; 
	// The point on the lens the ray passes through is in pLens
	Point2f pLens; 
	Float time;
};

class Camera
{

public:
	AnimatedTransform CameraToWorld;
	const Float shutterOpen, shutterClose;
	Film* film;
	const Medium* medium;


public:

	Camera(const AnimatedTransform& CameraToWorld, Float shutterOpen,
		Float shutterClose, Film* film, const Medium* medium);
	~Camera() {}

	virtual Float GenerateRay(const CameraSample& sample, Ray* ray) const = 0;

	// computes a main ray like GenerateRay() but
	// also computes the corresponding rays for pixels shifted one pixel in the x and y directions
	// on the film plane
	Float GenerateRayDifferential(const CameraSample& sample, RayDifferential* rd) const;


private:

};

