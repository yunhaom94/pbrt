
#include "core/filter.h"

Filter::Filter(const Vector2f& radius)
	: radius(radius),
	invRadius(Vector2f(1 / radius.x(), 1 / radius.y())) { }

Float BoxFilter::Evaluate(const Point2f& p) const
{
	return 1.0;
}

Float TriangleFilter::Evaluate(const Point2f& p) const
{
	return std::max((Float)0, radius.x() - std::abs(p.x())) *
		std::max((Float)0, radius.y() - std::abs(p.y()));
}

Float GaussianFilter::Evaluate(const Point2f& p) const
{
	return Gaussian(p.x(), expX) * Gaussian(p.y(), expY);
}

Float GaussianFilter::Gaussian(Float d, Float expv) const
{
	return std::max((Float)0, Float(std::exp(-alpha * d * d) - expv));
}

MitchellFilter::MitchellFilter(const Vector2f& radius, Float B, Float C) : 
	Filter(radius), B(B), C(C) {}

Float MitchellFilter::Evaluate(const Point2f& p) const 
{
	return Mitchell1D(p.x() * invRadius.x()) * Mitchell1D(p.y() * invRadius.y());
}

Float MitchellFilter::Mitchell1D(Float x) const 
{
	x = std::abs(2 * x);
	if (x > 1)
		return ((-B - 6 * C) * x * x * x + (6 * B + 30 * C) * x * x +
			(-12 * B - 48 * C) * x + (8 * B + 24 * C)) * (1.f / 6.f);
	else
		return ((12 - 9 * B - 6 * C) * x * x * x +
			(-18 + 12 * B + 6 * C) * x * x +
			(6 - 2 * B)) * (1.f / 6.f);
}

LanczosSincFilter::LanczosSincFilter(const Vector2f& radius, Float tau) : 
	Filter(radius), tau(tau) { }

Float LanczosSincFilter::Evaluate(const Point2f& p) const
{
	return WindowedSinc(p.x(), radius.x()) * WindowedSinc(p.y(), radius.y());
}
