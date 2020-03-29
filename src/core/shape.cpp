#include "core/shape.h"
#include "core/interaction.h"

Shape::Shape(const Transform* ObjectToWorld, const Transform* WorldToObject, bool reverseOrientation) : 
	ObjectToWorld(ObjectToWorld),
	WorldToObject(WorldToObject),
	reverseOrientation(reverseOrientation),
	transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {}

Bounds3d Shape::WorldBound() const {
	return (*ObjectToWorld)(ObjectBound());
}

bool Shape::IntersectP(const Ray& ray, bool testAlphaTexture) const
{
	double tHit = ray.tMax;
	SurfaceInteraction isect;
	return Intersect(ray, tHit, &isect, testAlphaTexture);
}
