#include "lights/distant.h"
#include "core/bounding_boxes.h"
#include "core/scene.h"

DistantLight::DistantLight(const Transform& LightToWorld,
	const Spectrum& L, const Vector3f& wLight)
	: Light((int)LightFlags::DeltaDirection, 
		LightToWorld,
		MediumInterface()),
	L(L), 
	wLight((LightToWorld(wLight)).normalized()) { }


void DistantLight::Preprocess(const Scene& scene) 
{
	scene.WorldBound().BoundingSphere(&worldCenter, &worldRadius);
}

Spectrum DistantLight::Sample_Li(const Interaction& ref,
	const Point2f& u, Vector3f* wi, Float* pdf,
	VisibilityTester* vis) const 
{
	*wi = wLight;
	*pdf = 1;
	Point3f pOutside = ref.p + wLight * (2 * worldRadius);
	*vis = VisibilityTester(ref, Interaction(pOutside, ref.time,
		mediumInterface));
	return L;
}

Spectrum DistantLight::Power() const 
{
	return L * Pi * worldRadius * worldRadius;
}