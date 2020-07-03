#pragma once

#include "core/pbrt.h"
#include "core/material.h"

// a defused material
class MatteMaterial : public Material
{
private:
	std::shared_ptr<Texture<Spectrum>> Kd;
	std::shared_ptr<Texture<Float>> sigma, bumpMap;

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

