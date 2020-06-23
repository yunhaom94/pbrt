#pragma once
#include "cameras/projective_camera.h"

class PerspectiveCamera : public ProjectiveCamera
{
public:
	
private:
	Vector3f dxCamera, dyCamera;

public:
	PerspectiveCamera(
		const AnimatedTransform& CameraToWorld,
		const Bounds2f& screenWindow, Float shutterOpen,
		Float shutterClose, Float lensRadius, Float focalDistance,
		Float fov, Film* film, const Medium* medium);

	Float GenerateRay(const CameraSample& sample, Ray* ray) const;

	Float GenerateRayDifferential(const CameraSample& sample,
		RayDifferential* ray) const;
};