#pragma once

#include "core/pbrt.h"
#include "core/memory.h"

// indicates whether the surface intersection was
// found along a path starting from the camera or 
// one starting from a light source
enum class TransportMode 
{
	Radiance, Importance 
};

class Material
{
public:
	// determining the reflective properties at the point and initializing the
	// BSDFs at si
	// Input: SurfaceInteraction* si
	virtual void ComputeScatteringFunctions(SurfaceInteraction* si,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const = 0;

	void Bump(const std::shared_ptr<Texture<Float>>& d, SurfaceInteraction* si) const;
};

