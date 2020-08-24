#include "filters/gaussian.h"

Float GaussianFilter::Evaluate(const Point2f& p) const
{
	return Gaussian(p.x(), expX) * Gaussian(p.y(), expY);
}

Float GaussianFilter::Gaussian(Float d, Float expv) const
{
	return std::max((Float)0, Float(std::exp(-alpha * d * d) - expv));
}