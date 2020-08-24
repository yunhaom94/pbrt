#pragma once

#include "core/pbrt.h"

class Filter 
{
public:
	// radius of filtering on x and y 
	// (actually a box from -x to + x and -y to +y)
	const Vector2f radius;
	const Vector2f invRadius;

public:
	Filter(const Vector2f& radius);
	
	// gives position of the sample point relative to the
	// center of the filter. Return the filter value
	virtual Float Evaluate(const Point2f& p) const = 0;

};