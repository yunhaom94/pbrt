#pragma once

#include "core/pbrt.h"
#include "core/primitive.h"
#include "core/bounding_boxes.h"

struct KdAccelNode 
{
	
	// anonymous union, same bits but different types to do things
	// save memory space and speed up
	// maybe change later TODO
	union 
	{
		Float split; // Interior
		int onePrimitive; // Leaf
		int primitiveIndicesOffset; // Leaf
	};
	union 
	{
		int flags; // Both
		int nPrims; // Leaf
		int aboveChild; // Interior
	};

	void InitLeaf(int* primNums, int np,
		std::vector<int>* primitiveIndices);

	void InitInterior(int axis, int ac, Float s);

	Float SplitPos() const { return split; }
	int nPrimitives() const { return nPrims >> 2; }
	int SplitAxis() const { return flags & 3; }
	bool IsLeaf() const { return (flags & 3) == 3; }
	int AboveChild() const { return aboveChild >> 2; }

};

enum class EdgeType { Start, End };
struct BoundEdge {
	// BoundEdge Public Methods
	BoundEdge() {}
	BoundEdge(Float t, int primNum, bool starting) : t(t), primNum(primNum) {
		type = starting ? EdgeType::Start : EdgeType::End;
	}
	Float t;
	int primNum;
	EdgeType type;
};

struct KdToDo {
	const KdAccelNode* node;
	Float tMin, tMax;
};

class KdTreeAccel : public Aggregate
{
public:
	const int isectCost, traversalCost, maxPrims;
	const Float emptyBonus;
	std::vector<std::shared_ptr<Primitive>> primitives;

private:
	std::vector<int> primitiveIndices;
	KdAccelNode* nodes;
	int nAllocedNodes, nextFreeNode;
	Bounds3f bounds;


public:
	KdTreeAccel(const std::vector<std::shared_ptr<Primitive>>& p,
		int isectCost, int traversalCost, Float emptyBonus,
		int maxPrims, int maxDepth);

	void buildTree(int nodeNum, const Bounds3f& nodeBounds,
		const std::vector<Bounds3f>& allPrimBounds, int* primNums,
		int nPrimitives, int depth,
		const std::unique_ptr<BoundEdge[]> edges[3],
		int* prims0, int* prims1, int badRefines=0);

	bool Intersect(const Ray& ray,
		SurfaceInteraction* isect) const;



};