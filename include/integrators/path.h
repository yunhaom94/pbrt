#pragma once

#include "integrators/sampler_integrator.h"

// integrator that does path tracing
class PathIntegrator : public SamplerIntegrator 
{
private:
	const int maxDepth;

public:
	PathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler);

	Spectrum Li(const RayDifferential& r, const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) const;
};

