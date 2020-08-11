#pragma once

#include "core/pbrt.h"
#include "core/camera.h"
#include "core/transform.h"


class ProjectiveCamera : public Camera
{

protected:
	Transform CameraToScreen, RasterToCamera;
	Transform ScreenToRaster, RasterToScreen;
	Float lensRadius, focalDistance;

public:
	ProjectiveCamera(const AnimatedTransform& CameraToWorld,
		const Transform& CameraToScreen,
		// screen space ratio, see below:
		/*    
		float frame = Float(film.fullResolution.x()) / Float(film.fullResolution.y());
		Bounds2f screen;
		if (frame > 1.f) {
			screen.pMin.x() = -frame;
			screen.pMax.x() = frame;
			screen.pMin.y() = -1.f;
			screen.pMax.y() = 1.f;
		}
		else {
			screen.pMin.x() = -1.f;
			screen.pMax.x() = 1.f;
			screen.pMin.y() = -1.f / frame;
			screen.pMax.y() = 1.f / frame;
		}
		*/
		const Bounds2f& screenWindow,
		Float shutterOpen, 
		Float shutterClose,
		Float lensr, 
		Float focald,
		Film* film, 
		const Medium* medium);
	

};

