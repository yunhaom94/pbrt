#include <math.h>

#include "shapes/sphere.h"
#include "core/transform.h"
#include "core/ray.h"
#include "core/interaction.h"
#include "utlis/utlis.h"

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
	//TODO: Ray ray = (*WorldToObject)(r, &oErr, &dErr);
	
	Ray ray = (*WorldToObject)(r);

	Eigen::Vector3d e = ray.o;
	Eigen::Vector3d d = ray.d;

	double t = 0, min_t = 0;
	double t1 = 0, t2 = 0;
	Eigen::Vector3d n;

	double discriminant = std::pow((d.dot(e)), 2) - (d.dot(d)) * ((e).dot(e) - std::pow(radius, 2));

	if (discriminant < 0)
		return false;

	// no intersection
	if (discriminant < 0)
		return false;
	// one intersection
	else if (discriminant == 0)
	{
		t = (0 - (d.dot(e))) / (d.dot(d));
		if (t < min_t)
			return false;
	}
	else if (discriminant > 0)
	{
		t1 = (0 - (d.dot(e)) + sqrt(discriminant)) / (d.dot(d));
		t2 = (0 - (d.dot(e)) - sqrt(discriminant)) / (d.dot(d));

		// if both behind screen
		if ((t1 < min_t && t2 < min_t) || (t1 > ray.tMax && t2 > ray.tMax) ||
			(t1 < min_t && t2 > ray.tMax) || (t2 < min_t && t1 > ray.tMax))
			return false;

		if (t1 < min_t)
			t = t2;
		else if (t2 < min_t)
			t = t1;
		if (t1 > ray.tMax)
			t = t2;
		else if (t2 > ray.tMax)
			t = t1;
		else if (t1 < t2)
			t = t1;
		else if (t1 >= t2)
			t = t2;
	}

	// compute hit position
	Eigen::Vector3d pHit = ray(t);

	// TODO: Refine sphere intersection point 225
	
	if (pHit.x() == 0 && pHit.y() == 0) 
		pHit[0] = 1e-5f * radius;
	
	double phi = atan2(pHit.y(), pHit.x());
	if (phi < 0) 
		phi += 2 * M_PI;

	// Test sphere intersection against clipping parameters
	// if zMIN/zMAX within radius and hit point within zMIN/MAX
	if ((zMin > -radius && pHit.z() < zMin) ||
		(zMax < radius && pHit.z() > zMax) || phi > phiMax)
	{
		if (t == t1)
			t = t2;
		else
			t = t1;

		pHit = ray(t);

		if (pHit.x() == 0 && pHit.y() == 0)
			pHit[0] = 1e-5f * radius;

		phi = atan2(pHit.y(), pHit.x());
		if (phi < 0)
			phi += 2 * M_PI;

		if ((zMin > -radius && pHit.z() < zMin) ||
			(zMax < radius && pHit.z() > zMax) || phi > phiMax)
			return false;
	}

	// get the 2D u,v coordinates
	// idea similiar to hw5
	double u = phi / phiMax;
	double theta = acos(Clamp(pHit.z() / radius, -1, 1));
	double v = (theta - thetaMin) / (thetaMax - thetaMin);

	//dpdu, dpdv calculation, see page 138
	double zRadius = std::sqrt(pHit.x() * pHit.x() + pHit.y() * pHit.y());
	double invZRadius = 1 / zRadius;
	double cosPhi = pHit.x() * invZRadius;
	double sinPhi = pHit.y() * invZRadius;
	Eigen::Vector3d dpdu(-phiMax * pHit.y(), phiMax * pHit.x(), 0);
	Eigen::Vector3d dpdv = (thetaMax - thetaMin) * Eigen::Vector3d(pHit.z() * cosPhi, pHit.z() * sinPhi, -radius * sin(theta));

	n = 2 * (e + t * d).normalized();
	isect->n = n;

	// TODO: partial derivaives of normal
	Eigen::Vector3d dndu(0, 0, 0);
	Eigen::Vector3d dndv(0, 0, 0);

	// TODO:
	Eigen::Vector3d pError(0, 0, 0);

	*isect = (*ObjectToWorld)(
		SurfaceInteraction(pHit, pError, Eigen::Vector2d(u, v), -ray.d, dpdu, dpdv,
			dndu, dndv, ray.time, this));

	tHit = t;

	return true;
}

double Sphere::Area() const
{
	return phiMax * radius * (zMax - zMin);
}
