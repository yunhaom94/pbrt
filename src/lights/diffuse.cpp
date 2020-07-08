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