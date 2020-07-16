#include "integrators/directlight.h"
#include "core/sampler.h"
#include "core/scene.h"
#include "core/sampling.h"

DirectLightingIntegrator::DirectLightingIntegrator(LightStrategy strategy, 
	int maxDepth,
	std::shared_ptr<const Camera> camera, 
	std::shared_ptr<Sampler> sampler) : 
	SamplerIntegrator(camera, sampler),
	strategy(strategy),
	maxDepth(maxDepth) {}

void DirectLightingIntegrator::Preprocess(const Scene& scene,
	Sampler& sampler)
{
	if (strategy == LightStrategy::UniformSampleAll)
	{
		for (const auto& light : scene.lights)
			nLightSamples.push_back(sampler.RoundCount(light->nSamples));

		for (int i = 0; i < maxDepth; ++i)
		{
			for (size_t j = 0; j < scene.lights.size(); ++j)
			{
				sampler.Request2DArray(nLightSamples[j]);
				sampler.Request2DArray(nLightSamples[j]);
			}
		}
	}
}

Spectrum DirectLightingIntegrator::Li(const RayDifferential& ray,
    const Scene& scene, Sampler& sampler,
    MemoryArena& arena, int depth) const
{
    Spectrum L(0.f);
    // Find closest ray intersection or return background radiance
    SurfaceInteraction isect;
    if (!scene.Intersect(ray, &isect)) {
        for (const auto& light : scene.lights) L += light->Le(ray);
        return L;
    }

    // Compute scattering functions for surface interaction
    isect.ComputeScatteringFunctions(ray, arena);
    if (!isect.bsdf)
        return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
    Vector3f wo = isect.wo;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);
    if (scene.lights.size() > 0) {
        // Compute direct lighting for _DirectLightingIntegrator_ integrator
        if (strategy == LightStrategy::UniformSampleAll)
            L += UniformSampleAllLights(isect, scene, arena, sampler,
                nLightSamples);
        else
            L += UniformSampleOneLight(isect, scene, arena, sampler);
    }
    if (depth + 1 < maxDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
        L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
    }
    return L;
}