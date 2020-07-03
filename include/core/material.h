#pragma once

#include "core/pbrt.h"
#include "core/bxdf.h"
#include "core/interaction.h"
#include "utlis/utlis.h"





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

