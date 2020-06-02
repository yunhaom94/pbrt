#include "core/scene.h"

Scene::Scene(std::shared_ptr<Primitive> aggregate, const std::vector<std::shared_ptr<Light>>& lights) : lights(lights), aggregate(aggregate)
{
	worldBound = *aggregate->WorldBound(); // TODO:
	for (const auto& light : lights)
		light->Preprocess(*this);
}

bool Scene::Intersect(const Ray& ray, SurfaceInteraction* isect) const {
	return aggregate->Intersect(ray, isect);
}

bool Scene::IntersectP(const Ray& ray) const
{
	return aggregate->IntersectP(ray);
}
