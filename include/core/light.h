#pragma once

#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/ray.h"
#include "core/medium_interface.h"
#include "core/interaction.h"
#include "core/transform.h"

enum class LightFlags : int 
{
	DeltaPosition = 1, DeltaDirection = 2, Area = 4, Infinite = 8
};

class Light
{
public:
	// light type
	const int flags; 
	// used by area light to compute soft shadow
	const int nSamples;
	const MediumInterface mediumInterface;

protected:
	const Transform LightToWorld, WorldToLight;

public:

	Light(int flags, const Transform& LightToWorld, const MediumInterface& mediumInterface, int nSamples = 1);

	virtual ~Light() {}
	
	virtual void Preprocess(const Scene& scene) {}
	
	Spectrum Le(Ray r) const { return Spectrum(0.0); }

	// returns the radiance arriving at the given point of interaction
	// VisibilityTester used to holds info on shadow ray
	// For area light:
	// u is the point on the surface of light
	// pdf output parameter stores the probability density for the light sample that was taken
	virtual Spectrum Sample_Li(const Interaction& ref, const Point2f& u,
		Vector3f* wi, Float* pdf, VisibilityTester* vis) const = 0;

	virtual Spectrum Power() const = 0;

	virtual Float Pdf_Li(const Interaction& ref,
		const Vector3f& wi) const = 0;



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
	// each end point of the shadow ray to be traced
	Interaction p0, p1;

public:
	VisibilityTester() {}

	VisibilityTester(const Interaction& p0, const Interaction& p1)
		: p0(p0), p1(p1) { }

	~VisibilityTester() {}

	const Interaction& P0() const { return p0; }
	const Interaction& P1() const { return p1; }

	// traces a shadow ray between 2 point and returns a Boolean result
	bool Unoccluded(const Scene& scene) const;
	
	// if ray pass through scattering medium
	// it accumulates the ray¡¯s transmittance
	Spectrum Tr(const Scene& scene,
		Sampler& sampler) const;

};

inline bool IsDeltaLight(int flags) {
	return flags & (int)LightFlags::DeltaPosition ||
		flags & (int)LightFlags::DeltaDirection;
}