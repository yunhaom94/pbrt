#include "core/shape.h"

Shape::Shape(const Transform* ObjectToWorld, const Transform* WorldToObject, bool reverseOrientation) : 
	ObjectToWorld(ObjectToWorld),
	WorldToObject(WorldToObject),
	reverseOrientation(reverseOrientation),
	transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {}

Bounds3d Shape::WorldBound() const {
	return (*ObjectToWorld)(ObjectBound());
}