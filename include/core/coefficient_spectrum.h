#pragma once
#include "core/pbrt.h"

template <int nSpectrumSamples> class CoefficientSpectrum
{
public:
	// template for static array allocation, i think
	double c[nSpectrumSamples];

public:
	CoefficientSpectrum(double v = 0) 
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] = v;
	}

	// sum of 2 spectrum are the sum of coefficient
	// same for other operations
	CoefficientSpectrum& operator+=(const CoefficientSpectrum& s2) 
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] += s2.c[i];
		return *this;
	}

	CoefficientSpectrum operator+(const CoefficientSpectrum& s2) const
	{
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] += s2.c[i];
		return ret;
	}

	CoefficientSpectrum& operator-=(const CoefficientSpectrum& s2) 
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] -= s2.c[i];
		return *this;
	}

	CoefficientSpectrum operator-(const CoefficientSpectrum& s2) const 
	{
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] -= s2.c[i];
		return ret;
	}

	CoefficientSpectrum& operator*=(const CoefficientSpectrum& s2)
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] *= s2.c[i];
		return *this;
	}

	CoefficientSpectrum operator*(const CoefficientSpectrum& s2) const 
	{
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] *= s2.c[i];
		return ret;
	}

	// TODO: division and unary negation

	friend CoefficientSpectrum Sqrt(const CoefficientSpectrum& s) {
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = sqrt(s.c[i]);
		return ret;
	}

	bool IsBlack() const
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] != 0.) return false;
		return true;

	}

	// linearly interpolate, why?
	inline Spectrum Lerp(double t, const Spectrum& s1, const Spectrum& s2)
	{
		return (1 - t) * s1 + t * s2;
	}

	// TODO: p316
	bool HasNaNs() const;

protected:




};

