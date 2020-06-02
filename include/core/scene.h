#pragma once

#include "core/pbrt.h"
#include "core/bounding_boxes.h"
#include "core/primitive.h"
#include "core/light.h"




class Scene
{

public:
	// list of light sources
	std::vector<std::shared_ptr<Light>> lights;

private:
	// kinda like the "top level" of a tree toward all geometries in 
	// the scene
	std::shared_ptr<Primitive> aggregate;
	Bounds3d worldBound;

public:

	Scene(std::shared_ptr<Primitive> aggregate, const std::vector<std::shared_ptr<Light>>& lights);
	~Scene() {}

	const Bounds3d& WorldBound() const { return worldBound; }

	// traces the given ray into the scene and returns a Boolean value indicating
	// if the ray intersected any of the primitives
	// out: insect - the closest intersection point
	bool Intersect(const Ray& ray, SurfaceInteraction* isect) const;

	// only check if there is interaction, used for shadow rays
	bool IntersectP(const Ray& ray) const;


private:

};

