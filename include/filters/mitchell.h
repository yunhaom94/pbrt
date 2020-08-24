#pragma once
#include "filters/box.h"

class MitchellFilter : public Filter
{
private:
	const Float B, C;

public:
	MitchellFilter(const Vector2f& radius, Float B, Float C);

	Float Evaluate(const Point2f& p) const;

	Float Mitchell1D(Float x) const;

};