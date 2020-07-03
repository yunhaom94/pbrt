#pragma once

#include "core/material.h"

class FourierMaterial : public Material
{
private:
	std::shared_ptr<Texture<Float>> bumpMap;
	FourierBSDFTable bsdfTable;

public:

	FourierMaterial(const std::string& filename, const std::shared_ptr<Texture<Float>>& bumpMap);

	~FourierMaterial() {}

	void ComputeScatteringFunctions(SurfaceInteraction* si, MemoryArena& arena, TransportMode mode, bool allowMultipleLobes) const;

};

