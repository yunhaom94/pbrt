#pragma once
#include "core/pbrt.h"

// This is the class with all the BSDF algorithms
// it is able to actually render a image
class Integrator
{
public:
	virtual void Render(const Scene& scene) = 0;

private:

};

Spectrum UniformSampleAllLights(const Interaction& it,
	const Scene& scene, MemoryArena& arena, Sampler& sampler,
	const std::vector<int>& nLightSamples, bool handleMedia = false);

Spectrum UniformSampleOneLight(const Interaction& it,
	const Scene& scene, MemoryArena& arena, Sampler& sampler,
	bool handleMedia = false);

Spectrum EstimateDirect(const Interaction& it,
	const Point2f& uScattering, const Light& light,
	const Point2f& uLight, const Scene& scene, Sampler& sampler,
	MemoryArena& arena, bool handleMedia = false, bool specular = false);