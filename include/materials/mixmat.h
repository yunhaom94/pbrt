#pragma once

#include "core/material.h"

class MixMaterial : public Material
{
private:
	std::shared_ptr<Material> m1, m2;
	std::shared_ptr<Texture<Spectrum>> scale;
	friend class MixMaterial;

public:
	MixMaterial(const std::shared_ptr<Material>& m1,
		const std::shared_ptr<Material>& m2,
		const std::shared_ptr<Texture<Spectrum>>& scale)
		: m1(m1), m2(m2), scale(scale) { }

	~MixMaterial() {}

	void ComputeScatteringFunctions(SurfaceInteraction* si,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const;

private:

};
