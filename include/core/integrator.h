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

	// do pre-computation on the scene after scene is initialized
	virtual void Preprocess(const Scene& scene, Sampler& sampler) { }

	/*
	* Rendering Loop:
	* 1. The Sampler provides a sequence of sample values, one for each image sample to be taken.
	* 2. The Camera turns a sample into a corresponding ray from the film plane
	* 3. Li() method implementation computes the radiance along that ray arriving at the film.
	*    (Li() == RayColor())
	* 4. The sample and its radiance are given to the Film, which stores their contribution in an image
	*/
	void Render(const Scene& scene);

	/* Get radiance of the camera ray, similar to Ray_Colour()
	* How it works:
	* 1. Use ray to intersect scene to get a surface interaction
	* 2. Compute material properties at interaction point in a BSDF
	* 3. Uses the light to determine the illumination
	* 4. Uses above to calculate the final radiance back to camera along the ray
	*/
	virtual Spectrum Li(const RayDifferential& ray, const Scene& scene,
		Sampler& sampler, MemoryArena& arena, int depth = 0) const = 0;

	Spectrum SpecularReflect(
		const RayDifferential& ray, const SurfaceInteraction& isect,
		const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) const;

	Spectrum SpecularTransmit(
		const RayDifferential& ray, const SurfaceInteraction& isect,
		const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) const;

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