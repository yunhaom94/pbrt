#pragma once
#include "core/light.h"

class SpotLight : public Light
{
private:
	const Point3f pLight;
	const Spectrum I;
	// angle where lights start to diminish
	// angle to where lights completely off
	const Float cosTotalWidth, cosFalloffStart;

public:
	SpotLight(const Transform& LightToWorld,
		const MediumInterface& mediumInterface, const Spectrum& I,
		Float totalWidth, Float falloffStart);
	Spectrum Sample_Li(const Interaction& ref, const Point2f& u, Vector3f* wi, Float* pdf, VisibilityTester* vis) const;
	Float Falloff(const Vector3f& w) const;
	Spectrum Power() const;
};