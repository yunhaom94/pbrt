#pragma once

#include "core/pbrt.h"
#include "core/integrator.h"

// its rendering process is driven by a stream of
// samples from a Sampler
class SamplerIntegrator : public Integrator
{
public:

private:
	// the sampler used to generate sampling rays
	std::shared_ptr<Sampler> sampler;

protected:
	std::shared_ptr<const Camera> camera;

public:
	SamplerIntegrator(std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler)
		: camera(camera), sampler(sampler) { }

	~SamplerIntegrator() {}

	// waiting to be implemented
	virtual void Preprocess(const Scene& scene, Sampler& sampler) { }

	// actually calculate light stuff along the ray,
	// similar to RayColour() in 418 hw
	// return the spectrograph of a ray of light
	virtual Spectrum Li(const RayDifferential& ray, const Scene& scene,
		Sampler& sampler, MemoryArena& arena, int depth = 0) const = 0;

	void Render(const Scene& scene);

	Spectrum SpecularReflect(
		const RayDifferential& ray, const SurfaceInteraction& isect,
		const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) const;

	Spectrum SpecularTransmit(
		const RayDifferential& ray, const SurfaceInteraction& isect,
		const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) const;

private:

};

