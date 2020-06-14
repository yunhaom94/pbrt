#include "core/shape.h"
#include "core/interaction.h"
#include "core/bounding_boxes.h"
#include "core/transform.h"



Shape::Shape(const Transform* ObjectToWorld, const Transform* WorldToObject, bool reverseOrientation) : 
	ObjectToWorld(ObjectToWorld),
	WorldToObject(WorldToObject),
	reverseOrientation(reverseOrientation),
	transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {}

Bounds3f Shape::WorldBound() const
{
	return (*ObjectToWorld)(ObjectBound());
}

bool Shape::IntersectP(const Ray& ray, bool testAlphaTexture) const
{
	Float tHit = ray.tMax;
	SurfaceInteraction isect;
	return Intersect(ray, &tHit, &isect, testAlphaTexture);
}

