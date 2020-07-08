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

};