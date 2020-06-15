#pragma once

#include "core/pbrt.h"

// a combination of shapes and materials
// the bridge between the geometry processing and shading
class Primitive
{
	// TODO:
public:
	Primitive() {}
	~Primitive() {}

	virtual Bounds3f WorldBound() const = 0;
	
	virtual bool Intersect(const Ray& r, SurfaceInteraction*) const = 0;
	virtual bool IntersectP(const Ray& r) const = 0;
	
	// GetAreaLight(), returns a pointer to the AreaLight that describes the
	// primitive¡¯s emission distribution, if the primitive is itself a light source.
	// otherwise nullptr
	virtual const AreaLight* GetAreaLight() const = 0;
	
	// returns a pointer to the material instance assigned to the primitive.
	virtual const Material* GetMaterial() const = 0;

	// init bsdf functions
	virtual void ComputeScatteringFunctions(SurfaceInteraction* isect,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const = 0;

private:

};

