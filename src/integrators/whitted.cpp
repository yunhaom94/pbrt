#include "integrators/whitted.h"
#include "core/spectrum.h"
#include "core/interaction.h"
#include "core/scene.h"
#include "core/ray.h"
#include "core/visibility_tester.h"
#include "core/sampler.h"
#include "core/memory.h"

Spectrum WhittedIntegrator::Li(const RayDifferential& ray,
	const Scene& scene, Sampler& sampler, MemoryArena& arena,
	int depth) const 
{

	// the total result light
	Spectrum L(0.0);

	// doing first hit
	SurfaceInteraction isect;
	if (!scene.Intersect(ray, &isect)) 
	{
		// if no hit, just get lights from light source
		for (const auto& light : scene.lights)
			L += light->Le(ray);
		return L;
	}

	// surface normal
	Normal3f n = isect.shading.n;
	// normalized vector toward eye
	Vector3f wo = isect.wo;

	// BSDF
	isect.ComputeScatteringFunctions(ray, arena);

	// get light if surface is a light source
	// 0 if not
	L += isect.Le(wo);

	// gather light from all light sources

	for (const auto& light : scene.lights) {
		// normalized vector toward light source
		Vector3f wi;

		// probability density function of light contribution if sample
		Float pdf;

		// use for shadow ray
		VisibilityTester visibility;

		// light density(in pdf) on wi direction
		Spectrum Li = light->Sample_Li(isect, sampler.Get2D(), &wi,
			&pdf, &visibility);

		// if there is no light coming
		if (Li.IsBlack() || pdf == 0) 
			continue;
		
		// do bsdf then multiply by consine (with dot product)
		// to get actual color
		// if light ray is blocked by some object, then skip (for shadow)
		Spectrum f = isect.bsdf->f(wo, wi);
		if (!f.IsBlack() && visibility.Unoccluded(scene))
			L += f * Li * std::abs(wi.dot(n)) * (1.0 / pdf);
	}

	// reflection and refraction
	if (depth + 1 < maxDepth)
	{
		L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
		L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
	}


	return L;
}