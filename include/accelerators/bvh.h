#pragma once

#include "core/pbrt.h"
#include "core/primitive.h"
#include "core/bounding_boxes.h"

// top to bottom: better tree/slower build to worse tree but better build
// SAH = surface area heuristic
// HLBVH = efficiently bvh
// Middle = middle
// EqualCounts = 418 A4
enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

// the box itself
struct BVHPrimitiveInfo
{
	size_t primitiveNumber;
	Bounds3f bounds;
	Point3f centroid;

	BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f& bounds)
		: primitiveNumber(primitiveNumber), bounds(bounds),
		centroid(0.5 * bounds.pMin + 0.5 * bounds.pMax) { }
	~BVHPrimitiveInfo() {}

};

struct BVHBuildNode 
{
	
	Bounds3f bounds;
	BVHBuildNode* children[2];
	int splitAxis, firstPrimOffset, nPrimitives;


	void InitLeaf(int first, int n, const Bounds3f& b) 
	{
		firstPrimOffset = first; // used for array representation of the tree
		nPrimitives = n; 
		bounds = b;
		children[0] = children[1] = nullptr;
	}

	void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1)
	{
		children[0] = c0;
		children[1] = c1;
		bounds = Union(c0->bounds, c1->bounds);
		splitAxis = axis;
		nPrimitives = 0;
	}
};

class BVHAccel
{
public:

	const int maxPrimsInNode;
	const SplitMethod splitMethod;
	std::vector<std::shared_ptr<Primitive>> primitives;


public:
	BVHAccel(const std::vector<std::shared_ptr<Primitive>>& p,
		int maxPrimsInNode, SplitMethod splitMethod);
	~BVHAccel();

	BVHBuildNode* recursiveBuild(MemoryArena& arena,
		std::vector<BVHPrimitiveInfo>& primitiveInfo, int start,
		int end, int* totalNodes,
		std::vector<std::shared_ptr<Primitive>>& orderedPrims);

	// Linear time build
	// TODO: read bout it p270
	BVHBuildNode* HLBVHBuild(MemoryArena& arena,
		const std::vector<BVHPrimitiveInfo>& primitiveInfo,
		int* totalNodes,
		std::vector<std::shared_ptr<Primitive>>& orderedPrims) const;

private:

};

