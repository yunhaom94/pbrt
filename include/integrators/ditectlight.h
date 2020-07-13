#pragma once
#include "integrators/sampler_integrator.h"

enum class LightStrategy { UniformSampleAll, UniformSampleOne };

class DirectLightingIntegrator : public SamplerIntegrator 
{
private:
	const LightStrategy strategy;
	const int maxDepth;
	std::vector<int> nLightSamples;

public:
    DirectLightingIntegrator(LightStrategy strategy, int maxDepth,
        std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler);

	void Preprocess(const Scene& scene, Sampler& sampler);
        
	Spectrum Li(const RayDifferential& ray, const Scene& scene,
		Sampler& sampler, MemoryArena& arena, int depth) const;

	
};