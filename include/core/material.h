#pragma once

#include "core/pbrt.h"
#include "core/memory.h"


enum class TransportMode 
{
	Radiance, Importance 
};

class Material
{
public:
	virtual void ComputeScatteringFunctions(SurfaceInteraction* si,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const = 0;

	void Bump(const std::shared_ptr<Texture<Float>>& d, SurfaceInteraction* si) const;
};

