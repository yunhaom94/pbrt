#pragma once

#include "core/pbrt.h"
#include "cameras/projective.h"
#include "core/transform.h"

class OrthographicCamera : public ProjectiveCamera 
{

private:
	Vector3f dxCamera, dyCamera;

public:
	OrthographicCamera(const AnimatedTransform& CameraToWorld,
		const Bounds2f& screenWindow, Float shutterOpen,
		Float shutterClose, Float lensRadius, Float focalDistance,
		Film* film, const Medium* medium);

	Float GenerateRay(const CameraSample& sample,
		Ray* ray) const;

	Float GenerateRayDifferential(
		const CameraSample& sample, RayDifferential* ray) const;

private:

};