#include "integrators/path.h"
#include "core/spectrum.h"
#include "core/ray.h"
#include "core/interaction.h"
#include "core/scene.h"
#include "core/reflection.h"
#include "core/sampler.h"

PathIntegrator::PathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
	std::shared_ptr<Sampler> sampler)
	: SamplerIntegrator(camera, sampler), maxDepth(maxDepth) { }

Spectrum PathIntegrator::Li(const RayDifferential& r, const Scene& scene,
	Sampler& sampler, MemoryArena& arena, int depth) const 
{
	Spectrum L(0.f), beta(1.f);
	RayDifferential ray(r);
	bool specularBounce = false;
	for (int bounces = 0; ; ++bounces)
	{
		// Intersect ray with sceneand store intersection in isect
		SurfaceInteraction isect;
		bool foundIntersection = scene.Intersect(ray, &isect);

		//Possibly add emitted light at intersection
		if (bounces == 0 || specularBounce) 
		{
			if (foundIntersection)
				L += beta * isect.Le(-ray.d);
			else
				for (const auto& light : scene.lights)
					L += beta * light->Le(ray);
		}

		// Terminate path if ray escaped or maxDepth was reached
		if (!foundIntersection || bounces >= maxDepth)
			break;

		// Compute scattering functions and skip over medium boundaries
		isect.ComputeScatteringFunctions(ray, arena, true);
		if (!isect.bsdf) {
			ray = isect.SpawnRay(ray.d);
			bounces--;
			continue;
		}
		
		// Sample illumination from lights to find path contribution 
		L += beta * UniformSampleOneLight(isect, scene, arena, sampler);

		// Sample BSDF to get new path direction 
		Vector3f wo = -ray.d, wi;
		Float pdf;
		BxDFType flags;
		Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(),
			&pdf, BSDF_ALL, &flags);
		if (f.IsBlack() || pdf == 0.f)
			break;
		beta *= f * std::abs(wi.dot(isect.shading.n)) / pdf;
		specularBounce = (flags & BSDF_SPECULAR) != 0;
		ray = isect.SpawnRay(wi);

		// TODO: Account for subsurface scattering, if applicable 915

		// Possibly terminate the path with Russian roulette 
		if (bounces > 3) {
			Float q = std::max((Float).05, 1 - beta.y());
			if (sampler.Get1D() < q)
				break;
			beta /= 1 - q;
		}

	}
	return L;
}