#pragma once

#include "core/light.h"


// directional light
class DistantLight : public Light 
{
private:
	const Spectrum L;
	const Vector3f wLight;
	Point3f worldCenter;
	Float worldRadius;

public:
	DistantLight(const Transform& LightToWorld, const Spectrum& L, const Vector3f& wLight);

	void Preprocess(const Scene& scene);

	Spectrum Sample_Li(const Interaction& ref, const Point2f& u, Vector3f* wi, Float* pdf, VisibilityTester* vis) const;

	Spectrum Power() const;

};