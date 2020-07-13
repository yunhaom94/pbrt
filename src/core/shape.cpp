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

Interaction Shape::Sample(const Interaction& ref, const Point2f& u) const
{
	return Sample(u);
}

Float Shape::Pdf(const Interaction&) const
{
	return 1.0 / Area();
}

Float Shape::Pdf(const Interaction& ref, const Vector3f& wi) const
{

	Ray ray = ref.SpawnRay(wi);
	Float tHit;
	SurfaceInteraction isectLight;
	if (!Intersect(ray, &tHit, &isectLight, false)) 
		return 0;

	Float pdf = (ref.p, isectLight.p).squaredNorm() /
		(std::abs(isectLight.n.dot(-wi)) * Area());

	return pdf;
}

