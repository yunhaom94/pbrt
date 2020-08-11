#pragma once
#include "cameras/perspective.h"
#include "core/ray.h"
#include "core/sampling.h"

PerspectiveCamera::PerspectiveCamera(
	const AnimatedTransform& CameraToWorld,
	const Bounds2f& screenWindow, Float shutterOpen,
	Float shutterClose, Float lensRadius, Float focalDistance,
	Float fov, Film* film, const Medium* medium)
	: ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.0),
		screenWindow, shutterOpen, shutterClose,
		lensRadius, focalDistance, film, medium)
{
	dxCamera = (RasterToCamera(Point3f(1, 0, 0)) -
		RasterToCamera(Point3f(0, 0, 0)));
	dyCamera = (RasterToCamera(Point3f(0, 1, 0)) -
		RasterToCamera(Point3f(0, 0, 0)));

	//TODO: Compute image plane bounds at z = 1 for PerspectiveCamera 951


}

inline Float PerspectiveCamera::GenerateRay(const CameraSample& sample, Ray* ray) const
{
	Point3f pFilm = Point3f(sample.pFilm.x(), sample.pFilm.y(), 0);
	Point3f pCamera = RasterToCamera(pFilm);

	// since all rays starts from origin in camera space, the direction
	// is the location of pixel of film once is transformed into camera space
	*ray = Ray(Point3f(0, 0, 0), Vector3f(pCamera).normalized());

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

inline Float PerspectiveCamera::GenerateRayDifferential(const CameraSample& sample, RayDifferential* ray) const
{
	Point3f pFilm = Point3f(sample.pFilm.x(), sample.pFilm.y(), 0);
	Point3f pCamera = RasterToCamera(pFilm);

	*ray = RayDifferential(Ray(pCamera, Vector3f(0, 0, 1)));

	if (lensRadius > 0) {
		Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
		Float ft = focalDistance / ray->d.z();
		Point3f pFocus = (*ray)(ft);
		ray->o = Point3f(pLens.x(), pLens.y(), 0);
		ray->d = (pFocus - ray->o).normalized();
	}
	else {
		ray->rxOrigin = ray->ryOrigin = ray->o;
		ray->rxDirection = (Vector3f(pCamera) + dxCamera).normalized();
		ray->ryDirection = (Vector3f(pCamera) + dyCamera).normalized();
	}

	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
	ray->hasDifferentials = true;
	ray->medium = medium;
	*ray = CameraToWorld(*ray);
	return 1;
}
