#pragma once

#include "core/pbrt.h"
#include "core/primitive.h"

// SAH = urface area heuristic
// HLBVH = 
// Middle = middle
// EqualCounts = 418 A4
enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

class BVHAccel
{
public:

	std::vector<std::shared_ptr<Primitive>> primitives;
	int maxPrimsInNode;
	SplitMethod splitMethod;


public:
	BVHAccel(const std::vector<std::shared_ptr<Primitive>>& p,
		int maxPrimsInNode, SplitMethod splitMethod);
	~BVHAccel() {}

private:

};

