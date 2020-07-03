#pragma once

#include "core/pbrt.h"
#include "core/bxdf.h"
#include "core/interaction.h"
#include "core/memory.h"
#include "core/texture.h"
#include "core/spectrum.h"
#include "utlis/utlis.h"


class Material
{
public:
	virtual void ComputeScatteringFunctions(SurfaceInteraction* si,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const = 0;

	void Bump(const std::shared_ptr<Texture<Float>>& d, SurfaceInteraction* si) const;
};

