#pragma once

#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/ray_differential.h"

// TODO: p741, move implementation and extra includes to .cpp
class Light
{
public:
	Light() {}
	~Light() {}
	
	void Preprocess(const Scene& scene) {}
	Spectrum Le(Ray r) { return Spectrum(0.0); }

	virtual Spectrum Sample_Li(const Interaction& ref, const Point2f& u,
		Vector3f* wi, Float* pdf,
		VisibilityTester* vis) const = 0;

private:

};

