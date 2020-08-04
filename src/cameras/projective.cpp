#include "cameras/projective.h"
#include "core/film.h"

ProjectiveCamera::ProjectiveCamera(const AnimatedTransform& CameraToWorld,
	const Transform& CameraToScreen, 
	const Bounds2f& screenWindow,
	Float shutterOpen,
	Float shutterClose, 
	Float lensr,
	Float focald,
	Film* film,
	const Medium* medium)
	: Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
	CameraToScreen(CameraToScreen)
{
	lensRadius = lensr;
	focalDistance = focald;

	ScreenToRaster =
		Scale(film->fullResolution.x(), film->fullResolution.y(), 1) *
		Scale(1 / (screenWindow.pMax.x() - screenWindow.pMin.x()),
			1 / (screenWindow.pMin.y() - screenWindow.pMax.y()), 1) *
		Translate(Vector3f(-screenWindow.pMin.x(), -screenWindow.pMax.y(), 0));
	RasterToScreen = Inverse(ScreenToRaster);
	RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;


}