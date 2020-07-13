#include "lights/diffuse.h"
#include "core/shape.h"

DiffuseAreaLight::DiffuseAreaLight(const Transform& LightToWorld,
	const MediumInterface& mediumInterface, const Spectrum& Lemit,
	int nSamples, const std::shared_ptr<Shape>& shape)
	: AreaLight(LightToWorld, mediumInterface, nSamples), Lemit(Lemit),
	shape(shape), area(shape->Area()) { }


Spectrum DiffuseAreaLight::L(const Interaction& intr, const Vector3f& w) const {
	return intr.n.dot(w) > 0.0 ? Lemit : Spectrum(0.0);
}

Spectrum DiffuseAreaLight::Power() const {
	return Lemit * area * Pi;
}

Spectrum DiffuseAreaLight::Sample_Li(const Interaction& ref, const Point2f& u, Vector3f* wi, Float* pdf, VisibilityTester* vis) const
{
	Interaction pShape = shape->Sample(ref, u);
	pShape.mediumInterface = mediumInterface;
	*wi = (pShape.p - ref.p).normalized();
	*pdf = shape->Pdf(ref, *wi);
	*vis = VisibilityTester(ref, pShape);
	return L(pShape, -*wi);
}

Float DiffuseAreaLight::Pdf_Li(const Interaction& ref, const Vector3f& wi) const
{
	return shape->Pdf(ref, wi);
	return Float();
}
