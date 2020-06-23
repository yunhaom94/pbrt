#include "core/camera.h"
#include "core/ray.h"

Camera::Camera(const AnimatedTransform& CameraToWorld, Float shutterOpen, Float shutterClose, Film* film, const Medium* medium) :
	CameraToWorld(CameraToWorld), shutterOpen(shutterOpen), shutterClose(shutterClose), film(film), medium(medium) {}

Float Camera::GenerateRayDifferential(const CameraSample& sample, RayDifferential* rd) const
{
	Float wt = GenerateRay(sample, rd);

	CameraSample sshift = sample;
	sshift.pFilm.x()++;
	Ray rx;
	Float wtx = GenerateRay(sshift, &rx);
	if (wtx == 0)
		return 0;
	rd->rxOrigin = rx.o;
	rd->rxDirection = rx.d;

	CameraSample yshift = sample;
	yshift.pFilm.y()++;
	Ray ry;
	Float wty = GenerateRay(yshift, &ry);
	if (wty == 0)
		return 0;
	rd->ryOrigin = ry.o;
	rd->ryDirection = ry.d;

	rd->hasDifferentials = true;
	return wt;
}
