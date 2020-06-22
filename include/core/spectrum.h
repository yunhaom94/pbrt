#pragma once
#include "core/pbrt.h"
#include "utlis/utlis.h"

// Spectrum Utility Declarations

// convert xyz to rgb
inline void XYZToRGB(const Float xyz[3], Float rgb[3]) 
{
	rgb[0] = 3.240479 * xyz[0] - 1.537150 * xyz[1] - 0.498535 * xyz[2];
	rgb[1] = -0.969256 * xyz[0] + 1.875991 * xyz[1] + 0.041556 * xyz[2];
	rgb[2] = 0.055648 * xyz[0] - 0.204043 * xyz[1] + 1.057311 * xyz[2];
}

inline void RGBToXYZ(const Float rgb[3], Float xyz[3]) 
{
	xyz[0] = 0.412453 * rgb[0] + 0.357580 * rgb[1] + 0.180423 * rgb[2];
	xyz[1] = 0.212671 * rgb[0] + 0.715160 * rgb[1] + 0.072169 * rgb[2];
	xyz[2] = 0.019334 * rgb[0] + 0.119193 * rgb[1] + 0.950227 * rgb[2];
}

template <int nSpectrumSamples> 
class CoefficientSpectrum
{
public:
	// template for static array allocation, i think
	Float c[nSpectrumSamples];
	static const int nSamples = nSpectrumSamples;

public:
	CoefficientSpectrum(Float v = 0)
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

	CoefficientSpectrum& operator/=(const CoefficientSpectrum& s2) const
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] /= s2.c[i];
		return *this;
	}

	CoefficientSpectrum operator/(const CoefficientSpectrum& s2) const
	{
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] /= s2.c[i];
		return ret;
	}

	bool operator==(const CoefficientSpectrum& s2) const
	{

		for (int i = 0; i < nSpectrumSamples; ++i)
		{
			if (c[i] != s2.c[i])
				return false;
		}

		return true
	}

	bool operator!=(const CoefficientSpectrum& s2) const
	{

		return (*this == s2);
	}

	CoefficientSpectrum operator*(Float a) const
	{
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] *= a;

		return ret;
	}
	CoefficientSpectrum& operator*=(Float a) 
	{
		for (int i = 0; i < nSpectrumSamples; ++i) 
			c[i] *= a;

		return *this;
	}

	// inline friend: https://stackoverflow.com/questions/381164/friend-and-inline-method-whats-the-point
	friend inline CoefficientSpectrum operator*(Float a,
		const CoefficientSpectrum& s)
	{
		return s * a;
	}

	friend CoefficientSpectrum Sqrt(const CoefficientSpectrum& s) {
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::sqrt(s.c[i]);
		return ret;
	}

	bool IsBlack() const
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
		{
			if (c[i] != 0)
				return false;
		}
		return true;

	}

	// linearly interpolate
	inline Spectrum Lerp(Float t, const Spectrum& s1, const Spectrum& s2)
	{
		return (1 - t) * s1 + t * s2;
	}

	CoefficientSpectrum Clamp(Float low = 0, Float high = Infinity) const
	{
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = ::Clamp(c[i], low, high);// global name-space
		return ret;
	}

	// check if anything is not a number
	// in case accident division by 0
	bool HasNaNs() const 
	{
		for (int i = 0; i < nSpectrumSamples; ++i)
		{
			if (std::isnan(c[i])) 
				return true;
		}
		return false;
	}

	Float& operator[](int i)
	{
		return c[i];
	}

protected:


};

// sample spectrum definitions:
static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60;

static const int nCIESamples = 471;
static const Float CIE_Y_integral = 106.856895;
// these should be given
extern const Float CIE_X[nCIESamples];
extern const Float CIE_Y[nCIESamples];
extern const Float CIE_Z[nCIESamples];
extern const Float CIE_lambda[nCIESamples];

static const int nRGB2SpectSamples = 32;
extern const Float RGB2SpectLambda[nRGB2SpectSamples];
extern const Float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const Float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const Float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const Float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const Float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const Float RGBRefl2SpectBlue[nRGB2SpectSamples];

extern const Float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const Float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const Float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const Float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const Float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const Float RGBIllum2SpectBlue[nRGB2SpectSamples];

enum class SpectrumType { Reflectance, Illuminant };

// helper functions
bool SpectrumSamplesSorted(const Float* lambda, const Float* v, int n);
void SortSpectrumSamples(Float* lambda, Float* v, int n);
Float AverageSpectrumSamples(const Float* lambda, const Float* vals, int n, Float lambdaStart, Float lambdaEnd);
Float InterpolateSpectrumSamples(const Float* lambda, const Float* vals,
	int n, Float l);


class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples>
{
public:

private:
	// X,Y,Z response curves
	static SampledSpectrum X, Y, Z;

	// used to convert RGB back to SPD p328
	static SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
	static SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
	static SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
	static SampledSpectrum rgbRefl2SpectBlue;

	static SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
	static SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
	static SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
	static SampledSpectrum rgbIllum2SpectBlue;

public:
	// given a provided spectral data as a set of (¦Ëi , vi)
	// where the ith sample has some value vi at wavelength ¦Ë
	// convert to SampledSpectrum
	static SampledSpectrum FromSampled(const Float* lambda, const Float* v, int n)
	{

		if (!SpectrumSamplesSorted(lambda, v, n)) {
			// create a new copy to this data
			std::vector<Float> slambda(&lambda[0], &lambda[n]);
			std::vector<Float> sv(&v[0], &v[n]);
			// sort it
			SortSpectrumSamples(&slambda[0], &sv[0], n);
			return FromSampled(&slambda[0], &sv[0], n);
		}

		SampledSpectrum r;
		// spiting the spectrum into n equal segments,
		// each starts with lambda 0 and ends with lambda 1
		// then get average value of that interval from input
		for (int i = 0; i < nSpectralSamples; ++i) {
			Float lambda0 = ::Lerp(Float(i) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			Float lambda1 = ::Lerp(Float(i + 1) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
		}

		return r;
	}


	// compute XYZ matching curves for this SampleSpectrum
	static void Init() 
	{
		for (int i = 0; i < nSpectralSamples; ++i) 
		{
			Float wl0 = ::Lerp(Float(i) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			Float wl1 = ::Lerp(Float(i + 1) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples,
				wl0, wl1);
			Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples,
				wl0, wl1);
			Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples,
				wl0, wl1);
		}

		// Compute RGB to spectrum functions for _SampledSpectrum_
		for (int i = 0; i < nSpectralSamples; ++i) 
		{
			Float wl0 = ::Lerp(Float(i) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			Float wl1 = ::Lerp(Float(i + 1) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			rgbRefl2SpectWhite.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
					nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectCyan.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
					nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectMagenta.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
					nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectYellow.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
					nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectRed.c[i] = AverageSpectrumSamples(
				RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectGreen.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
					nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectBlue.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
					nRGB2SpectSamples, wl0, wl1);

			rgbIllum2SpectWhite.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
					nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectCyan.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
					nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectMagenta.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
					nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectYellow.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
					nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectRed.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
					nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectGreen.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
					nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectBlue.c[i] =
				AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
					nRGB2SpectSamples, wl0, wl1);
		}
	}

	static SampledSpectrum FromRGB(const Float rgb[3], SpectrumType type);

	static SampledSpectrum FromXYZ(const Float xyz[3],
		SpectrumType type = SpectrumType::Reflectance) 
	{
		Float rgb[3];
		XYZToRGB(xyz, rgb);
		return FromRGB(rgb, type);
	}


	SampledSpectrum(Float v = 0) : CoefficientSpectrum(v) { }

	SampledSpectrum(const CoefficientSpectrum<nSpectralSamples>& v)
		: CoefficientSpectrum<nSpectralSamples>(v) {}

	SampledSpectrum(const RGBSpectrum& r, SpectrumType t);

	// convert spectrum to XYZ color 
	// used to convert to RGB color later
	void ToXYZ(Float xyz[3]) const
	{
		xyz[0] = xyz[1] = xyz[2] = 0.f;
		for (int i = 0; i < nSpectralSamples; ++i)
		{
			xyz[0] += X.c[i] * c[i];
			xyz[1] += Y.c[i] * c[i];
			xyz[2] += Z.c[i] * c[i];
		}
		Float scale = Float(sampledLambdaEnd - sampledLambdaStart) /
			Float(CIE_Y_integral * nSpectralSamples);
		xyz[0] *= scale;
		xyz[1] *= scale;
		xyz[2] *= scale;

	}

	// The y coefficient of XYZ color is important
	Float y() const 
	{
		Float yy = 0.f;
		for (int i = 0; i < nSpectralSamples; ++i)
			yy += Y.c[i] * c[i];
		
		return yy * Float(sampledLambdaEnd - sampledLambdaStart) /
			Float(CIE_Y_integral * nSpectralSamples);
	}

	void ToRGB(Float rgb[3]) const
	{
		Float xyz[3];
		ToXYZ(xyz);
		XYZToRGB(xyz, rgb);
	}


	

	RGBSpectrum ToRGBSpectrum() const;

private:
};


class RGBSpectrum : public CoefficientSpectrum<3>
{
public:
	RGBSpectrum(Float v = 0) : CoefficientSpectrum<3>(v) { }
	RGBSpectrum(const CoefficientSpectrum<3>& v) : CoefficientSpectrum<3>(v) { }

	// converting RGB values to spectrum is trivial
	// same signature for compile time polymorphism
	static RGBSpectrum FromRGB(const Float rgb[3], SpectrumType type = SpectrumType::Reflectance)
	{
		RGBSpectrum s;
		s.c[0] = rgb[0];
		s.c[1] = rgb[1];
		s.c[2] = rgb[2];
		return s;
	}

	// same as above
	void ToRGB(Float* rgb) const
	{
		rgb[0] = c[0];
		rgb[1] = c[1];
		rgb[2] = c[2];
	}

	const RGBSpectrum& ToRGBSpectrum() const

	{
		return *this;
	}

	void ToXYZ(Float xyz[3]) const 
	{ 
		RGBToXYZ(c, xyz); 
	}

	static RGBSpectrum FromXYZ(const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) 
	{
		RGBSpectrum r;
		XYZToRGB(xyz, r.c);
		return r;
	}

	Float y() const 
	{
		const Float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
		return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
	}

	// converts the spectrum to XYZand then to RGB.
	static RGBSpectrum FromSampled(const Float* lambda, const Float* v,
		int n)
	{

		if (!SpectrumSamplesSorted(lambda, v, n)) 
		{
			// create a new copy to this data
			std::vector<Float> slambda(&lambda[0], &lambda[n]);
			std::vector<Float> sv(&v[0], &v[n]);
			// sort it
			SortSpectrumSamples(&slambda[0], &sv[0], n);
			return FromSampled(&slambda[0], &sv[0], n);
		}

		Float xyz[3] = { 0, 0, 0 };
		for (int i = 0; i < nCIESamples; ++i)
		{
			Float val = InterpolateSpectrumSamples(lambda, v, n, CIE_lambda[i]);
			xyz[0] += val * CIE_X[i];
			xyz[1] += val * CIE_Y[i];
			xyz[2] += val * CIE_Z[i];
		}
		Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
			Float(CIE_Y_integral * nCIESamples);
		xyz[0] *= scale;
		xyz[1] *= scale;
		xyz[2] *= scale;
		return FromXYZ(xyz);
	}


};