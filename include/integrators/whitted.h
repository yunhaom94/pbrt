#pragma once

#include "core/pbrt.h"
#include "core/integrator.h"

// a simple integrator without direct lighting
class WhittedIntegrator : public SamplerIntegrator
{

private:
	// max recursive depth for reflection
	const int maxDepth;

public:
	WhittedIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler)
		: SamplerIntegrator(camera, sampler), maxDepth(maxDepth) { }

	Spectrum Li(const RayDifferential& ray, const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) const;

	~WhittedIntegrator() {}

private:

	

};

