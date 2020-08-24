#pragma once
#include "filters/box.h"

// weight on a Gaussian function
class GaussianFilter : public Filter
{
private:
	const Float alpha;
	const Float expX, expY;

public:
	GaussianFilter(const Vector2f& radius, Float alpha)
		: Filter(radius), alpha(alpha),
		expX(std::exp(-alpha * radius.x() * radius.x())),
		expY(std::exp(-alpha * radius.y() * radius.y())) { }

	Float Evaluate(const Point2f& p) const;

	Float Gaussian(Float d, Float expv) const;

};