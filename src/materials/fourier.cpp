#include "materials/fourier.h"

FourierMaterial::FourierMaterial(const std::string& filename,
	const std::shared_ptr<Texture<Float>>& bumpMap)
	: bumpMap(bumpMap) {
	FourierBSDFTable::Read(filename, &bsdfTable);
}

void FourierMaterial::ComputeScatteringFunctions(SurfaceInteraction* si,
	MemoryArena& arena, TransportMode mode,
	bool allowMultipleLobes) const
{

	if (bumpMap)
		Bump(bumpMap, si);
	si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
	si->bsdf->Add(ARENA_ALLOC(arena, FourierBSDF)(bsdfTable, mode));
}