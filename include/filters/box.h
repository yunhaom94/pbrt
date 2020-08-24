#pragma once
#include "core/filter.h"

// Equally weights all samples within a square region of the image.
class BoxFilter : public Filter
{
public:
	BoxFilter(const Vector2f& radius) : Filter(radius) { }
	Float Evaluate(const Point2f& p) const;
};