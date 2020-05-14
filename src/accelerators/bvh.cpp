#include "accelerators/bvh.h"

BVHAccel::BVHAccel(const std::vector<std::shared_ptr<Primitive>>& p,
	int maxPrimsInNode, SplitMethod splitMethod)
	: primitives(p), maxPrimsInNode(std::min(255, maxPrimsInNode)),
	splitMethod(splitMethod) {
	if (primitives.size() == 0)
		return;
	
}