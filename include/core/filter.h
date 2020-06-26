#pragma once

#include "core/pbrt.h"

class Filter 
{
public:
	const Vector2f radius, invRadius;

public:
	Filter(const Vector2f& radius);
	virtual Float Evaluate(const Point2f& p) const = 0;

};

class BoxFilter : public Filter 
{
public:
	BoxFilter(const Vector2f& radius) : Filter(radius) { }
	Float Evaluate(const Point2f& p) const;
};

class TriangleFilter : public Filter
{
public:
	TriangleFilter(const Vector2f& radius) : Filter(radius) { }
	Float Evaluate(const Point2f& p) const;
};

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

class MitchellFilter : public Filter 
{
private:
	const Float B, C;

public:
	MitchellFilter(const Vector2f& radius, Float B, Float C);

	Float Evaluate(const Point2f& p) const;

	Float Mitchell1D(Float x) const;

};

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