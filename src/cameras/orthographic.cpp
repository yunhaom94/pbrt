#include "cameras/orthographic.h"
#include "core/film.h"
#include "core/ray.h"
#include "core/sampling.h"

OrthographicCamera::OrthographicCamera(const AnimatedTransform& CameraToWorld,
	const Bounds2f& screenWindow, Float shutterOpen,
	Float shutterClose, Float lensRadius, Float focalDistance,
	Film* film, const Medium* medium)
	: ProjectiveCamera(CameraToWorld, Orthographic(0, 1),
		screenWindow, shutterOpen, shutterClose,
		lensRadius, focalDistance, film, medium)
{
	dxCamera = RasterToCamera(Vector3f(1, 0, 0));
	dyCamera = RasterToCamera(Vector3f(0, 1, 0));
}

Float OrthographicCamera::GenerateRay(const CameraSample& sample, Ray* ray) const
{
	Point3f pFilm = Point3f(sample.pFilm.x(), sample.pFilm.y(), 0);
	Point3f pCamera = RasterToCamera(pFilm);
	*ray = Ray(pCamera, Vector3f(0, 0, 1));

	if (lensRadius > 0)
	{
		Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
		Float ft = focalDistance / ray->d.z();
		Point3f pFocus = (*ray)(ft);
		ray->o = Point3f(pLens.x(), pLens.y(), 0);
		ray->d = (pFocus - ray->o).normalized();
	}

	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
	ray->medium = medium;
	*ray = CameraToWorld(*ray);
	return 1;
}

Float OrthographicCamera::GenerateRayDifferential(const CameraSample& sample, RayDifferential* ray) const
{
	Point3f pFilm = Point3f(sample.pFilm.x(), sample.pFilm.y(), 0);
	Point3f pCamera = RasterToCamera(pFilm);

	*ray = RayDifferential(Ray(pCamera, Vector3f(0, 0, 1)));

	if (lensRadius > 0) 
	{
		Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
		Float ft = focalDistance / ray->d.z();
		Point3f pFocus = (*ray)(ft);
		ray->o = Point3f(pLens.x(), pLens.y(), 0);
		ray->d = (pFocus - ray->o).normalized();
	}
	else 
	{
		ray->rxOrigin = ray->o + dxCamera;
		ray->ryOrigin = ray->o + dyCamera;
		ray->rxDirection = ray->ryDirection = ray->d;
	}

	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
	ray->hasDifferentials = true;
	ray->medium = medium;
	*ray = CameraToWorld(*ray);
	return 1;
}
