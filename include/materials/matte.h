#pragma once

#include "core/pbrt.h"
#include "core/material.h"

// a defuse only material
class MatteMaterial : public Material
{
private:
	// spectral diffuse reflection value
	std::shared_ptr<Texture<Spectrum>> Kd;
	// scalar roughness value 
	// (if sigma = 0, then lambertian otherwise orennayar)
	std::shared_ptr<Texture<Float>> sigma;
	std::shared_ptr<Texture<Float>> bumpMap;

public:
	MatteMaterial(const std::shared_ptr<Texture<Spectrum>>& Kd,
		const std::shared_ptr<Texture<Float>>& sigma,
		const std::shared_ptr<Texture<Float>>& bumpMap)
		: Kd(Kd), sigma(sigma), bumpMap(bumpMap) { }

	~MatteMaterial() {}

	void ComputeScatteringFunctions(SurfaceInteraction* si,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const;
};

//MatteMaterial* CreateMatteMaterial(const TextureParams& mp);