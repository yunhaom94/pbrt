#pragma once

#include "core/material.h"

// a mix of diffuse and specular
class PlasticMaterial : public Material 
{
private:
	std::shared_ptr<Texture<Spectrum>> Kd, Ks;
	std::shared_ptr<Texture<Float>> roughness, bumpMap;
	const bool remapRoughness;

public:
	PlasticMaterial(const std::shared_ptr<Texture<Spectrum>>& Kd,
		const std::shared_ptr<Texture<Spectrum>>& Ks,
		const std::shared_ptr<Texture<Float>>& roughness,
		const std::shared_ptr<Texture<Float>>& bumpMap,
		bool remapRoughness)
		: Kd(Kd), Ks(Ks), roughness(roughness), bumpMap(bumpMap),
		remapRoughness(remapRoughness) { }

	void ComputeScatteringFunctions(SurfaceInteraction* si, MemoryArena& arena, TransportMode mode, bool allowMultipleLobes) const;

};