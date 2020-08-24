#pragma once
#include "filters/box.h"

// samples at the filter center have a weight of 1, 
// and the weight falls off linearly to the square extent of the filter.
class TriangleFilter : public Filter
{
public:
	TriangleFilter(const Vector2f& radius) : Filter(radius) { }
	Float Evaluate(const Point2f& p) const;
};
