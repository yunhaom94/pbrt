#include "..\..\include\core\primitive.h"
#include "..\..\include\core\primitive.h"
#include "..\..\include\core\primitive.h"
#include "..\..\include\core\primitive.h"
#include "core/primitive.h"

GeometricPrimitive::GeometricPrimitive(std::shared_ptr<Shape> shape, std::shared_ptr<Material> material,
	std::shared_ptr<AreaLight> areaLight, MediumInterface mediumInterface) :
	shape(shape), material(material), areaLight(areaLight), mediumInterface(mediumInterface) {}

bool GeometricPrimitive::Intersect(const Ray& r,
	SurfaceInteraction* isect) const 
{
	Float tHit;
	if (!shape->Intersect(r, &tHit, isect))
		return false;
	r.tMax = tHit;
	isect->primitive = this;
	//TODO: Initialize SurfaceInteraction::mediumInterface after Shape intersection 685 
	return true;
}

Bounds3f GeometricPrimitive::WorldBound() const 
{ 
	return shape->WorldBound(); 
}

bool GeometricPrimitive::IntersectP(const Ray& r) const 
{
	return shape->IntersectP(r);
}

const AreaLight* GeometricPrimitive::GetAreaLight() const
{
	return areaLight.get();
}

const Material* GeometricPrimitive::GetMaterial() const
{
	return material.get();
}

void GeometricPrimitive::ComputeScatteringFunctions(SurfaceInteraction* isect,
	MemoryArena& arena, TransportMode mode, bool allowMultipleLobes) const
{
	if (material)
		material->ComputeScatteringFunctions(isect, arena, mode, allowMultipleLobes);
}



bool TransformedPrimitive::Intersect(const Ray& r, SurfaceInteraction* isect) const
{
	// transform ray to primitive's space
	Transform InterpolatedPrimToWorld;
	// some weird shit
	PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
	Ray ray = Inverse(InterpolatedPrimToWorld)(r);

	if (!primitive->Intersect(ray, isect))
		return false;
	r.tMax = ray.tMax;

	if (!InterpolatedPrimToWorld.IsIdentity())
		*isect = InterpolatedPrimToWorld(*isect);

	return true;
}

bool TransformedPrimitive::IntersectP(const Ray& r) const
{
	Transform InterpolatedPrimToWorld;

	PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
	Ray ray = Inverse(InterpolatedPrimToWorld)(r);

	if (!primitive->IntersectP(ray, isect))
		return false;

	return true;
}

AreaLight* TransformedPrimitive::GetAreaLight() const
{
	return nullptr;
}

Material* TransformedPrimitive::GetMaterial() const
{
	return nullptr;
}

void TransformedPrimitive::ComputeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allowMultipleLobes) const
{
	return;
}
