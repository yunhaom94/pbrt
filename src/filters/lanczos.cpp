#include "filters/lanczos.h"

LanczosSincFilter::LanczosSincFilter(const Vector2f& radius, Float tau) :
	Filter(radius), tau(tau) { }

Float LanczosSincFilter::Evaluate(const Point2f& p) const
{
	return WindowedSinc(p.x(), radius.x()) * WindowedSinc(p.y(), radius.y());
}
