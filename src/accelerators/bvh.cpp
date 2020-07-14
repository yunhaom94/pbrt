#include "accelerators/bvh.h"
#include "core/memory.h"

BVHAccel::BVHAccel(const std::vector<std::shared_ptr<Primitive>>& p,
	int maxPrimsInNode, SplitMethod splitMethod)
	: primitives(p), maxPrimsInNode(std::min(255, maxPrimsInNode)),
	splitMethod(splitMethod)
{
	if (primitives.size() == 0)
		return;

	// building bvh tree
	// the array form a the bvh
	std::vector<BVHPrimitiveInfo> primitiveInfo;
	primitiveInfo.reserve(primitives.size());

	for (size_t i = 0; i < primitives.size(); ++i)
		primitiveInfo[i] = { i, primitives[i]->WorldBound() };

	MemoryArena arena(1024 * 1024);
	int totalNodes = 0;
	// primitives shall be sorted, and the new pointer will be assigned to old
	// primitives pointer
	std::vector<std::shared_ptr<Primitive>> orderedPrims;
	BVHBuildNode* root;
	if (splitMethod == SplitMethod::HLBVH)
		root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
	else
		root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
			&totalNodes, orderedPrims);
	primitives.swap(orderedPrims);

}


BVHBuildNode* BVHAccel::recursiveBuild(MemoryArena& arena, std::vector<BVHPrimitiveInfo>& primitiveInfo,
	int start, int end, int* totalNodes,
	std::vector<std::shared_ptr<Primitive>>& orderedPrims)
{
	BVHBuildNode* node = arena.Alloc<BVHBuildNode>();
	(*totalNodes)++;

	// bound of all primitives
	Bounds3f bounds;
	for (int i = start; i < end; ++i)
		bounds = Union(bounds, primitiveInfo[i].bounds);

	int nPrimitives = end - start;
	// base case, create leaf
	if (nPrimitives == 1)
	{
		int firstPrimOffset = orderedPrims.size();
		for (int i = start; i < end; ++i)
		{
			int primNum = primitiveInfo[i].primitiveNumber;
			orderedPrims.push_back(primitives[primNum]);
		}
		node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
		return node;
	}
	else
	{
		// find the bound for all center points
		Bounds3f centroidBounds;
		for (int i = start; i < end; ++i)
			centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
		int dim = centroidBounds.MaximumExtent();
		int mid = (start + end) / 2;
		// just create a leaf
		if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim])
		{
			int firstPrimOffset = orderedPrims.size();
			for (int i = start; i < end; ++i) {
				int primNum = primitiveInfo[i].primitiveNumber;
				orderedPrims.push_back(primitives[primNum]);
			}
			node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
			return node;
		}
		else
		{
			switch (splitMethod)
			{
			case SplitMethod::Middle:
			{
				Float pmid = (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
				BVHPrimitiveInfo* midPtr =
					std::partition(&primitiveInfo[start], &primitiveInfo[end - 1] + 1,
						[dim, pmid](const BVHPrimitiveInfo& pi) {
					return pi.centroid[dim] < pmid;
				});
				mid = midPtr - &primitiveInfo[0];
				if (mid != start && mid != end)
					break;
			}
			case SplitMethod::EqualCounts:
			{
				mid = (start + end) / 2;
				std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
					&primitiveInfo[end - 1] + 1,
					[dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
					return a.centroid[dim] < b.centroid[dim];
				});
				break;
			}
			case SplitMethod::SAH:
			{
				if (nPrimitives <= 4)
				{	// equal division for small groups
					std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
						&primitiveInfo[end - 1] + 1,
						[dim](const BVHPrimitiveInfo& a,
							const BVHPrimitiveInfo& b) {
						return a.centroid[dim] < b.centroid[dim];
					});
					break;
				}
				else
				{
					// using something similar to a bucket sort?
					constexpr int nBuckets = 12;
					struct BucketInfo
					{
						int count = 0;
						Bounds3f bounds;
					};
					BucketInfo buckets[nBuckets];

					for (int i = start; i < end; ++i)
					{
						int b = nBuckets *
							centroidBounds.Offset(primitiveInfo[i].centroid)[dim];
						if (b == nBuckets) b = nBuckets - 1;
						buckets[b].count++;
						buckets[b].bounds = Union(buckets[b].bounds, primitiveInfo[i].bounds);
					}

					// compute cost for splitting
					Float cost[nBuckets - 1];
					for (int i = 0; i < nBuckets - 1; ++i)
					{
						Bounds3f b0, b1;
						int count0 = 0, count1 = 0;
						for (int j = 0; j <= i; ++j)
						{
							b0 = Union(b0, buckets[j].bounds);
							count0 += buckets[j].count;
						}
						for (int j = i + 1; j < nBuckets; ++j)
						{
							b1 = Union(b1, buckets[j].bounds);
							count1 += buckets[j].count;
						}
						cost[i] = 0.125 + (count0 * b0.SurfaceArea() +
							count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
					}

					Float minCost = cost[0];
					int minCostSplitBucket = 0;
					for (int i = 1; i < nBuckets - 1; ++i)
					{
						if (cost[i] < minCost)
						{
							minCost = cost[i];
							minCostSplitBucket = i;
						}
					}

					Float leafCost = nPrimitives;
					if (nPrimitives > maxPrimsInNode || minCost < leafCost)
					{
						BVHPrimitiveInfo* pmid = std::partition(&primitiveInfo[start],
							&primitiveInfo[end - 1] + 1,
							[=](const BVHPrimitiveInfo& pi) {
							int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
							if (b == nBuckets) b = nBuckets - 1;
							return b <= minCostSplitBucket;
						});
						mid = pmid - &primitiveInfo[0];
					}
					else
					{
						int firstPrimOffset = orderedPrims.size();
						for (int i = start; i < end; ++i) {
							int primNum = primitiveInfo[i].primitiveNumber;
							orderedPrims.push_back(primitives[primNum]);
						}
						node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
						return node;
					}

				}


				break;
			}
			}

			node->InitInterior(dim,
				recursiveBuild(arena, primitiveInfo, start, mid,
					totalNodes, orderedPrims),
				recursiveBuild(arena, primitiveInfo, mid, end,
					totalNodes, orderedPrims));
		}

	}
	return node;
}

BVHBuildNode* BVHAccel::HLBVHBuild(MemoryArena& arena, const std::vector<BVHPrimitiveInfo>& primitiveInfo, 
	int* totalNodes, std::vector<std::shared_ptr<Primitive>>& orderedPrims) const
{
	// TODO: p268
	return nullptr;
}
