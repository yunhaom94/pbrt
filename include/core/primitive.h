#pragma once

#include "core/pbrt.h"

class Primitive
{
	// TODO:
public:
	Primitive() {}
	~Primitive() {}
	void virtual IntersectP() = 0;
	Bounds3d* WorldBound() { return NULL; }
	bool Intersect(const Ray& ray, SurfaceInteraction* isect) { return false; }
	bool IntersectP(const Ray& ray) { return false; }

private:

};

