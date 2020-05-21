#pragma once

#include "cameras/camera.h"
#include "core/transform.h"

class ProjectiveCamera : public Camera
{

protected:
	Transform CameraToScreen, RasterToCamera;

public:
	ProjectiveCamera(const AnimatedTransform& CameraToWorld,
		const Transform& CameraToScreen, const Eigen::Vector3d& screenWindow,
		double shutterOpen, double shutterClose, double lensr, double focald,
		Film* film, const Medium* medium)
		: Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
		CameraToScreen(CameraToScreen);
	
protected:

};