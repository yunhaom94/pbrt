#pragma once

#include "core/pbrt.h"
#include "core/bounding_boxes.h"

class Scene
{

public:
	std::vector<std::shared_ptr<Light>> lights;
	std::shared_ptr<Primitive> aggregate;

private:
	Bounds3d worldBound;

public:
	// TODO: p24
	Scene(std::shared_ptr<Primitive> aggregate, const std::vector<std::shared_ptr<Light>>& lights) : lights(lights), aggregate(aggregate) {}
	~Scene();

	const Bounds3d& WorldBound() const { return worldBound; }
	void Intersect();

private:

};

