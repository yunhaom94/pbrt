#include "shapes/sphere.h"
#include "core/transform.h"


Sphere::Sphere(const Transform* ObjectToWorld, const Transform* WorldToObject, bool reverseOrientation, double radius, double zMin, double zMax, double phiMax)
	: Shape(ObjectToWorld, WorldToObject, reverseOrientation),
	radius(radius), 
	zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
	zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
	thetaMin(acos(Clamp(zMin / radius, -1, 1))),
	thetaMax(acos(Clamp(zMax / radius, -1, 1))),
	phiMax((Clamp(phiMax, 0, 2 * M_PI))) 
{
}


Bounds3d Sphere::ObjectBound() const
{
	return Bounds3d(Eigen::Vector3d(-radius, -radius, zMin),
		Eigen::Vector3d(radius, radius, zMax));
}

bool Sphere::Intersect(const Ray& r, double& tHit, SurfaceInteraction* isect, bool testAlphaTexture) const
{
	// offset in translation
	Eigen::Vector3d oErr, dErr;
	//Ray ray = (*WorldToObject)(r, &oErr, &dErr);


	return false;
}

double Sphere::Area() const
{
	return 0.0;
}
