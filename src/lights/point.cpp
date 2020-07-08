#include "lights/point.h"

PointLight::PointLight(const Transform& LightToWorld,
	const MediumInterface& mediumInterface,
	const Spectrum& I) :
	Light((int)LightFlags::DeltaPosition, LightToWorld, mediumInterface),
	pLight(LightToWorld(Point3f(0, 0, 0))), I(I) { }

Spectrum PointLight::Sample_Li(const Interaction& ref,
	const Point2f& u, Vector3f* wi, Float* pdf,
	VisibilityTester* vis) const
{
	*wi = (pLight - ref.p).normalized();
	*pdf = 1.f;
	*vis = VisibilityTester(ref, Interaction(pLight, ref.time,
		mediumInterface));
	return I / (pLight, ref.p).squaredNorm();
}

Spectrum PointLight::Power() const 
{
	return 4 * Pi * I;
}