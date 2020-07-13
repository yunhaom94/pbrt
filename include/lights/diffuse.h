#include "core/light.h"

class DiffuseAreaLight : public AreaLight
{
protected:
	const Spectrum Lemit;
	// this light emits from outward face of an shape 
	std::shared_ptr<Shape> shape;
	const Float area;

public:
	DiffuseAreaLight(const Transform& LightToWorld, const MediumInterface& mediumInterface, const Spectrum& Lemit, int nSamples, const std::shared_ptr<Shape>& shape);

	Spectrum L(const Interaction& intr, const Vector3f& w) const;

	Spectrum Power() const;

	Spectrum Sample_Li(const Interaction& ref,
		const Point2f& u, Vector3f* wi, Float* pdf,
		VisibilityTester* vis) const;

	Float Pdf_Li(const Interaction& ref,
		const Vector3f& wi) const;
};