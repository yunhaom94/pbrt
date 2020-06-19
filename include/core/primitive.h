#pragma once

#include "core/pbrt.h"
#include "core/medium_interface.h"
#include "core/transform.h"

// a combination of shapes and materials
// the bridge between the geometry processing and shading
class Primitive
{
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

// a single shape in scene
class GeometricPrimitive : public Primitive
{

public:
	std::shared_ptr<Shape> shape;
	std::shared_ptr<Material> material;
	// when its a light source
	std::shared_ptr<AreaLight> areaLight;
	MediumInterface mediumInterface;

public:
	GeometricPrimitive(std::shared_ptr<Shape> shape, std::shared_ptr<Material> material,
		std::shared_ptr<AreaLight> areaLight, MediumInterface mediumInterface);

	bool Intersect(const Ray& r, SurfaceInteraction* isect) const;
	bool IntersectP(const Ray& r) const;
	
	Bounds3f WorldBound() const;

	AreaLight* GetAreaLight() const;
	Material* GetMaterial() const;
	
	void ComputeScatteringFunctions(SurfaceInteraction* isect,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const;

private:

};

// used for object instancing, i.e same object, multiple places
// also used for animation
// takes a reference to the Primitive that represents the model,
// and the transformation that places it in the scene.
class TransformedPrimitive : public Primitive 
{

public:
	std::shared_ptr<Primitive> primitive;
	const AnimatedTransform PrimitiveToWorld;

public:
	TransformedPrimitive(std::shared_ptr<Primitive>& primitive,
		const AnimatedTransform& PrimitiveToWorld)
		: primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) { }

	bool Intersect(const Ray& r, SurfaceInteraction* isect) const;
	bool IntersectP(const Ray& r) const;

	// should not be called:
	AreaLight* GetAreaLight() const;
	Material* GetMaterial() const;
	void ComputeScatteringFunctions(SurfaceInteraction* isect,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const;
};

// the class for accelerate structures
// basically a interface to multiple primitives
class Aggregate : public Primitive 
{
public:


	// should not be called:
	AreaLight* GetAreaLight() const;
	Material* GetMaterial() const;
	void ComputeScatteringFunctions(SurfaceInteraction* isect,
		MemoryArena& arena, TransportMode mode,
		bool allowMultipleLobes) const;
};