#pragma once

#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/ray.h"
#include "core/medium_interface.h"
#include "core/interaction.h"
#include "core/transform.h"

enum class LightFlags : int {
	DeltaPosition = 1, DeltaDirection = 2, Area = 4, Infinite = 8
};

class Light
{
public:
	const int flags; // light type
	const int nSamples;
	const MediumInterface mediumInterface;

protected:
	const Transform LightToWorld, WorldToLight;

public:

	Light(int flags, const Transform& LightToWorld, const MediumInterface& mediumInterface, int nSamples = 1);

	virtual ~Light() {}
	
	virtual void Preprocess(const Scene& scene) {}
	
	Spectrum Le(Ray r) { return Spectrum(0.0); }

	virtual Spectrum Sample_Li(const Interaction& ref, const Point2f& u,
		Vector3f* wi, Float* pdf, VisibilityTester* vis) const = 0;

	virtual Spectrum Power() const = 0;

private:

};

class AreaLight : public Light
{
public:
	virtual Spectrum L(const Interaction& intr, const Vector3f& w) const = 0;

	AreaLight(const Transform& LightToWorld, const MediumInterface& medium, int nSamples = 1);

};

// Used for shadow ray
class VisibilityTester 
{
private:
	Interaction p0, p1;

public:
	VisibilityTester() {}

	VisibilityTester(const Interaction& p0, const Interaction& p1)
		: p0(p0), p1(p1) { }

	~VisibilityTester() {}

	const Interaction& P0() const { return p0; }
	const Interaction& P1() const { return p1; }

	bool Unoccluded(const Scene& scene) const;
	
	// if ray pass through scattering medium
	Spectrum Tr(const Scene& scene,
		Sampler& sampler) const;

};