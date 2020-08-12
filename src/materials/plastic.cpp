#include "materials/plastic.h"
#include "core/microfacet.h"
#include "core/interaction.h"
#include "core/texture.h"
#include "core/spectrum.h"
#include "core/reflection.h"
#include "utlis/utlis.h"

void PlasticMaterial::ComputeScatteringFunctions(
	SurfaceInteraction* si, MemoryArena& arena, TransportMode mode,
	bool allowMultipleLobes) const
{
	if (bumpMap)
		Bump(bumpMap, si);

	si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
	Spectrum kd = Kd->Evaluate(*si).Clamp();
	if (!kd.IsBlack())
		si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(kd));

	Spectrum ks = Ks->Evaluate(*si).Clamp();
	if (!ks.IsBlack())
	{
		Fresnel* fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, 1.5f);
		Float rough = roughness->Evaluate(*si);
		if (remapRoughness)
			rough = TrowbridgeReitzDistribution::RoughnessToAlpha(rough);
		
		MicrofacetDistribution* distrib = ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(rough, rough);
		BxDF* spec = ARENA_ALLOC(arena, MicrofacetReflection)(ks, distrib, fresnel);
		si->bsdf->Add(spec);
	}



}