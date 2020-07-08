#include "core/light.h"
#include "core/transform.h"
#include "core/scene.h"
#include "core/ray.h"
#include "core/medium.h"


Light::Light(int flags, 
	const Transform& LightToWorld,
	const MediumInterface& mediumInterface,
	int nSamples = 1)
	: flags(flags),
	nSamples(std::max(1, nSamples)),
	mediumInterface(mediumInterface),
	LightToWorld(LightToWorld),
	WorldToLight(Inverse(LightToWorld)) 
{
	// TODO: Warn if light has transformation with non-uniform scale	
}

AreaLight::AreaLight(const Transform& LightToWorld, const MediumInterface& medium,
	int nSamples)
	: Light((int)LightFlags::Area, LightToWorld, medium, nSamples) { }


bool VisibilityTester::Unoccluded(const Scene& scene) const {
	return !scene.IntersectP(p0.SpawnRayTo(p1));
}

Spectrum VisibilityTester::Tr(const Scene& scene, Sampler& sampler) const
{
	Ray ray(p0.SpawnRayTo(p1));
	Spectrum Tr(1.f);
	while (true) 
	{
		SurfaceInteraction isect;
		bool hitSurface = scene.Intersect(ray, &isect);
		if (hitSurface && isect.primitive->GetMaterial() != nullptr)
			return Spectrum(0.0);
		if (ray.medium)
			Tr *= ray.medium->Tr(ray, sampler);
		if (!hitSurface)
			break;
		ray = isect.SpawnRayTo(p1);
	}
	return Tr;
	return Spectrum();
}


