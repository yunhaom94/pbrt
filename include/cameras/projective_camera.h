#pragma once

#include "core/pbrt.h"
#include "core/camera.h"
#include "core/transform.h"


class ProjectiveCamera : public Camera
{

protected:
	Transform CameraToScreen, RasterToCamera;
	Float lensRadius, focalDistance;

public:
	ProjectiveCamera(const AnimatedTransform& CameraToWorld,
		const Transform& CameraToScreen,
		const Bounds2f& screenWindow,
		Float shutterOpen, 
		Float shutterClose,
		Float lensr, 
		Float focald,
		Film* film, 
		const Medium* medium);
	

};

