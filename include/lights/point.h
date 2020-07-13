#pragma once
#include "core/light.h"

class PointLight : public Light 
{
private:
	const Point3f pLight;
	const Spectrum I;

public:
	PointLight(const Transform& LightToWorld, const MediumInterface& mediumInterface, const Spectrum& I);

	Spectrum Sample_Li(const Interaction& ref, const Point2f& u, Vector3f* wi, Float* pdf, VisibilityTester* vis) const;

	Spectrum Power() const;

	Float Pdf_Li(const Interaction&, const Vector3f&) const;

	
};