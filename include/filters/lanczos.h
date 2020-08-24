#pragma once
#include "filters/box.h"

// a modified sinc filter, by multiplying
// a function that goes to 0 to get a finite extent
class LanczosSincFilter : public Filter
{
private:
	const Float tau;

public:
	LanczosSincFilter(const Vector2f& radius, Float tau);

	Float Evaluate(const Point2f& p) const;

	Float Sinc(Float x) const {
		x = std::abs(x);
		if (x < 1e-5) return 1;
		return std::sin(Pi * x) / (Pi * x);
	}

	Float WindowedSinc(Float x, Float radius) const {
		x = std::abs(x);
		if (x > radius) return 0;
		Float lanczos = Sinc(x / tau);
		return Sinc(x) * lanczos;
	}
};