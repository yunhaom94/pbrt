#pragma once

#include "core/camera.h"


// TODO: p376
class EnvironmentCamera : public Camera 
{
public:
	EnvironmentCamera(const AnimatedTransform& CameraToWorld,
		Float shutterOpen, Float shutterClose, Film* film,
		const Medium* medium)
		: Camera(CameraToWorld, shutterOpen, shutterClose, film, medium) {}

};