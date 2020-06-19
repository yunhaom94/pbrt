#pragma once

#include "core/pbrt.h"

class Material
{
public:
	Material() {}
	~Material() {}

	void ComputeScatteringFunctions(SurfaceInteraction* isect,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const {}

private:

};

