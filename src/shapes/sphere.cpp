#include "shapes/sphere.h"
#include "core/transform.h"
#include "core/ray.h"
#include "core/interaction.h"
#include "core/bounding_boxes.h"
#include "core/efloat.h"
#include "core/sampling.h"
#include "utlis/utlis.h"



Sphere::Sphere(const Transform* ObjectToWorld,
	const Transform* WorldToObject,
	bool reverseOrientation,
	Float radius,
	Float zMin,
	Float zMax,
	Float phiMax)
	: Shape(ObjectToWorld, WorldToObject, reverseOrientation),
	radius(radius), 
	zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
	zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
	thetaMin(std::acos(Clamp(zMin / radius, -1, 1))),
	thetaMax(std::acos(Clamp(zMax / radius, -1, 1))),
	phiMax(Radians(Clamp(phiMax, 0, 360)))
{
}


Bounds3f Sphere::ObjectBound() const
{
	return Bounds3f(Point3f(-radius, -radius, zMin),
		Point3f(radius, radius, zMax));
}

bool Sphere::Intersect(const Ray& r, Float* tHit,
	SurfaceInteraction* isect, bool testAlphaTexture) const
{


	// we do intersection in object space
	Vector3f oErr, dErr;
	Ray ray = (*WorldToObject)(r, &oErr, &dErr);

	// compute quadratic equation and check for solution:
	EFloat ox(ray.o.x(), oErr.x()), oy(ray.o.y(), oErr.y()), oz(ray.o.z(), oErr.z());
	EFloat dx(ray.d.x(), dErr.x()), dy(ray.d.y(), dErr.y()), dz(ray.d.z(), dErr.z());
	EFloat a = dx * dx + dy * dy + dz * dz;
	EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
	EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

	EFloat t0, t1;
	if (!Quadratic(a, b, c, &t0, &t1))
		return false;

	// since t0 <= t1
	if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0)
		return false;

	// check which one to pick:
	EFloat tShapeHit = t0;
	
	if (t0.LowerBound() <= 0) 
	{
		tShapeHit = t1;
		if (t1.UpperBound() > ray.tMax)
			return false;
	}


	// check for clipped sphere:
	Float phi;
	Point3f pHit = ray((Float)tShapeHit);

	// Refine sphere intersection point 225
	// finding hit point on the sphere coordinate
	if (pHit.x() == 0 && pHit.y() == 0)
		pHit.x() = 1e-5f * radius;

	phi = std::atan2(pHit.y(), pHit.x());

	if (phi < 0) 
		phi += 2 * Pi;

	if ((zMin > -radius && pHit.z() < zMin) ||
		(zMax < radius && pHit.z() > zMax) || phi > phiMax) 
	{
		if (tShapeHit == t1) 
			return false;
		if (t1.UpperBound() > ray.tMax) 
			return false;

		// do again for t1
		tShapeHit = t1;
		pHit = ray((Float)tShapeHit);

		// Refine sphere intersection point 225
		// finding hit point on the sphere coordinate
		if (pHit.x() == 0 && pHit.y() == 0)
			pHit.x() = 1e-5f * radius;

		if (phi < 0)
			phi += 2 * Pi;
		
		if ((zMin > -radius && pHit.z() < zMin) ||
			(zMax < radius && pHit.z() > zMax) || phi > phiMax)
			return false;
	}

	// calculate partial derivatives
	Float u = phi / phiMax;
	Float theta = std::acos(Clamp(pHit.z() / radius, -1, 1));
	Float v = (theta - thetaMin) / (thetaMax - thetaMin);

	Float zRadius = std::sqrt(pHit.x()* pHit.x()+ pHit.y() * pHit.y());
	Float invZRadius = 1.0 / zRadius;
	Float cosPhi = pHit.x()* invZRadius;
	Float sinPhi = pHit.y() * invZRadius;
	Vector3f dpdu(-phiMax * pHit.y(), phiMax * pHit.x(), 0);
	Vector3f dpdv = (thetaMax - thetaMin) *
		Vector3f(pHit.z() * cosPhi, pHit.z() * sinPhi, -radius * std::sin(theta));

	Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x(), pHit.y(), 0);
	Vector3f d2Pduv = (thetaMax - thetaMin) * pHit.z() * phiMax *
		Vector3f(-sinPhi, cosPhi, 0.);
	Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
		Vector3f(pHit.x(), pHit.y(), pHit.z());

	// PARTIAL DERIVATIVES OF NORMAL VECTORS
	Float E = dpdu.dot(dpdu);
	Float F = dpdu.dot(dpdv);
	Float G = dpdv.dot(dpdv);
	Vector3f N = (dpdu.cross(dpdv)).normalized();
	Float e = N.dot(d2Pduu);
	Float f = N.dot(d2Pduv);
	Float g = N.dot(d2Pdvv);

	Float invEGF2 = 1 / (E * G - F * F);
	Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
		(e * F - f * E) * invEGF2 * dpdv);
	Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
		(f * F - g * E) * invEGF2 * dpdv);

	// error shit
	Vector3f pError = gamma(5) * ((Vector3f)pHit.cwiseAbs());

	*isect = (*ObjectToWorld)(
		SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d, dpdu, dpdv,
			dndu, dndv, ray.time, this));

	*tHit = (Float)tShapeHit;

	return true;
}

bool Sphere::IntersectP(const Ray& r, bool testAlphaTexture) const
{
	// we do intersection in object space
	Vector3f oErr, dErr;
	Ray ray = (*WorldToObject)(r, &oErr, &dErr);

	// compute quadratic equation and check for solution:
	EFloat ox(ray.o.x(), oErr.x()), oy(ray.o.y(), oErr.y()), oz(ray.o.z(), oErr.z());
	EFloat dx(ray.d.x(), dErr.x()), dy(ray.d.y(), dErr.y()), dz(ray.d.z(), dErr.z());
	EFloat a = dx * dx + dy * dy + dz * dz;
	EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
	EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

	EFloat t0, t1;
	if (!Quadratic(a, b, c, &t0, &t1))
		return false;

	// since t0 <= t1
	if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0)
		return false;

	// check which one to pick:
	EFloat tShapeHit = t0;

	if (t0.LowerBound() <= 0)
	{
		tShapeHit = t1;
		if (t1.UpperBound() > ray.tMax)
			return false;
	}


	// check for clipped sphere:
	Float phi;
	Point3f pHit = ray((Float)tShapeHit);

	// TODO: Refine sphere intersection point 225
	// finding hit point on the sphere coordinate
	if (pHit.x() == 0 && pHit.y() == 0)
		pHit.x() = 1e-5f * radius;

	phi = std::atan2(pHit.y(), pHit.x());

	if (phi < 0)
		phi += 2 * Pi;

	if ((zMin > -radius && pHit.z() < zMin) ||
		(zMax < radius && pHit.z() > zMax) || phi > phiMax)
	{
		if (tShapeHit == t1)
			return false;
		if (t1.UpperBound() > ray.tMax)
			return false;

		// do again for t1
		tShapeHit = t1;
		pHit = ray((Float)tShapeHit);

		// TODO: Refine sphere intersection point 225
		// finding hit point on the sphere coordinate
		if (pHit.x() == 0 && pHit.y() == 0)
			pHit.x() = 1e-5f * radius;

		if (phi < 0)
			phi += 2 * Pi;

		if ((zMin > -radius && pHit.z() < zMin) ||
			(zMax < radius && pHit.z() > zMax) || phi > phiMax)
			return false;
	}


	return true;
}

double Sphere::Area() const
{
	return phiMax * radius * (zMax - zMin);
}

Interaction Sphere::Sample(const Point2f& u) const
{
	Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
	Interaction it;
	it.n = (*ObjectToWorld)(Normal3f(pObj.x(), pObj.y(), pObj.z()));
	it.n.normalize();
	
	if (reverseOrientation) 
		it.n *= -1;
	
	pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
	Vector3f pObjError = gamma(5) * ((Vector3f)pObj).cwiseAbs();
	
	it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
	return it;
}

Interaction Sphere::Sample(const Interaction& ref,
	const Point2f& u) const
{
	// Compute coordinate system for sphere sampling
	Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
	Vector3f wc = (pCenter - ref.p).normalized();
	Vector3f wcX, wcY;
	CoordinateSystem(wc, &wcX, &wcY);

	// Sample uniformly on sphere if p is inside it
	Point3f pOrigin = OffsetRayOrigin(ref.p, ref.pError, ref.n,
		pCenter - ref.p);
	if ((pOrigin - pCenter).squaredNorm() <= radius * radius)
		return Sample(u);

	// Sample sphere uniformly inside subtended cone

	Float sinThetaMax2 = radius * radius / (ref.p - pCenter).squaredNorm();
	Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
	Float cosTheta = (1 - u[0]) + u[0] * cosThetaMax;
	Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
	Float phi = u[1] * 2 * Pi;

	Float dc = Distance(ref.p, pCenter);
	Float ds = dc * cosTheta -
		std::sqrt(std::max((Float)0,
			radius * radius - dc * dc * sinTheta * sinTheta));
	Float cosAlpha = (dc * dc + radius * radius - ds * ds) /
		(2 * dc * radius);
	Float sinAlpha = std::sqrt(std::max((Float)0, 1 - cosAlpha * cosAlpha));

	Vector3f nObj = SphericalDirection(sinAlpha, cosAlpha, phi,
		-wcX, -wcY, -wc);
	Point3f pObj = radius * Point3f(nObj.x(), nObj.y(), nObj.z());

	Interaction it;
	pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
	Vector3f pObjError = gamma(5) * ((Vector3f)pObj).cwiseAbs();
	it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
	it.n = (*ObjectToWorld)(Normal3f(nObj));
	if (reverseOrientation) it.n *= -1;
	return it;
}

Float Sphere::Pdf(const Interaction& ref, const Vector3f& wi) const
{
	Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
	// Return uniform PDF if point is inside sphere 844

	Point3f pOrigin = OffsetRayOrigin(ref.p, ref.pError, ref.n,
		pCenter - ref.p);

	if ((pOrigin - pCenter).squaredNorm() <= radius * radius)
		return Shape::Pdf(ref, wi);

	Float sinThetaMax2 = radius * radius / (ref.p - pCenter).squaredNorm();
	Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
	return UniformConePdf(cosThetaMax);
}