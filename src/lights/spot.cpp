#include "lights/spot.h"
#include "..\..\include\lights\spot.h"

SpotLight::SpotLight(const Transform& LightToWorld,
	const MediumInterface& mediumInterface, const Spectrum& I,
	Float totalWidth, Float falloffStart)
	: Light((int)LightFlags::DeltaPosition, LightToWorld,
		mediumInterface),
	pLight(LightToWorld(Point3f(0, 0, 0))), I(I),
	cosTotalWidth(std::cos(Radians(totalWidth))),
	cosFalloffStart(std::cos(Radians(falloffStart))) { }

Spectrum SpotLight::Sample_Li(const Interaction& ref,
	const Point2f& u, Vector3f* wi, Float* pdf,
	VisibilityTester* vis) const {
	*wi = (pLight - ref.p).normalized();
	*pdf = 1.0;
	*vis = VisibilityTester(ref, Interaction(pLight, ref.time,
		mediumInterface));
	return I * Falloff(-*wi) / (pLight - ref.p).squaredNorm();
}

Float SpotLight::Falloff(const Vector3f& w) const 
{
	Vector3f wl = (WorldToLight(w)).normalized();
	Float cosTheta = wl.z();

	if (cosTheta < cosTotalWidth) 
		return 0;
	if (cosTheta > cosFalloffStart) 
		return 1;

	Float delta = (cosTheta - cosTotalWidth) /
		(cosFalloffStart - cosTotalWidth);
	return (delta * delta) * (delta * delta);
}

Spectrum SpotLight::Power() const 
{
	return I * 2 * Pi * (1 - 0.5 * (cosFalloffStart + cosTotalWidth));
}

