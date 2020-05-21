#pragma once

#include "core/coefficient_spectrum.h"
#include "utlis/utlis.h"

enum class SpectrumType { Reflectance, Illuminant };

static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60;

static const int nCIESamples = 471;
extern const double CIE_X[nCIESamples];
extern const double CIE_Y[nCIESamples];
extern const double CIE_Z[nCIESamples];
extern const double CIE_lambda[nCIESamples];
static const double CIE_Y_integral = 106.856895;


static bool SpectrumSamplesSorted(const double* lambda, const double* v, int n)
{
	double last_l = lambda[0];

	for (int i = 0; i < n; i++)
	{
		double l = lambda[i];
		// check if wavelength is not sorted
		if (l < last_l)
			return false;

		last_l = l;
	}

	return true;
}

static void SortSpectrumSamples(double* lambda, double* v, int n)
{

	// insertion sort 
	// TODO: change this if very slow

	int i, key, j;
	for (i = 1; i < n; i++)
	{
		key = lambda[i];
		j = i - 1;


		while (j >= 0 && lambda[j] > key)
		{
			lambda[j + 1] = lambda[j];
			v[j + 1] = v[j];
			j = j - 1;
		}
		lambda[j + 1] = key;
	}
}

// doing some average over area
static double AverageSpectrumSamples(const double* lambda, const double* vals, int n, double lambdaStart, double lambdaEnd)
{
	//boundary check
	if (lambdaEnd <= lambda[0]) return vals[0];
	if (lambdaStart >= lambda[n - 1]) return vals[n - 1];
	if (n == 1) return vals[0];

	double sum = 0;

	// Add contributions of constant segments before/after samples
	if (lambdaStart < lambda[0])
		sum += vals[0] * (lambda[0] - lambdaStart);
	if (lambdaEnd > lambda[n - 1])
		sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);

	int i = 0;
	while (lambdaStart > lambda[i + 1])
		++i;

	auto interp = [lambda, vals](double w, int i) {
		return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]),
			vals[i], vals[i + 1]);
	};
	for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) {
		double segLambdaStart = std::max(lambdaStart, lambda[i]);
		double segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
		sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
			(segLambdaEnd - segLambdaStart);
	}

	return sum / (lambdaEnd - lambdaStart);
}



static double InterpolateSpectrumSamples(const double* lambda, const double* vals,
	int n, double l) {
	if (l <= lambda[0]) return vals[0];
	if (l >= lambda[n - 1]) return vals[n - 1];
	int offset = FindInterval(n,
		[&](int index) { return lambda[index] <= l; });
	double t = (l - lambda[offset]) / (lambda[offset + 1] - lambda[offset]);
	return Lerp(t, vals[offset], vals[offset + 1]);
}


class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples>
{
public:
	SampledSpectrum(double v = 0) : CoefficientSpectrum(v) { }

	// 1. give a bunch of sample points
	// 2. connect them by lines
	// 3. calculate sample range and n samples split location
	// 4. calculate area of each sample segment
	// 5. average?
	static SampledSpectrum FromSampled(const double* lambda,
		const double* v, int n)
	{

		if (!SpectrumSamplesSorted(lambda, v, n)) {
			// create a new copy to this data
			std::vector<double> slambda(&lambda[0], &lambda[n]);
			std::vector<double> sv(&v[0], &v[n]);
			// sort it
			SortSpectrumSamples(&slambda[0], &sv[0], n);
			return FromSampled(&slambda[0], &sv[0], n);
		}

		SampledSpectrum r;
		// spiting the spectrum into n equal segments,
		// each starts with lambda 0 and ends with lambda 1
		for (int i = 0; i < nSpectralSamples; ++i) {
			double lambda0 = ::Lerp(double(i) / double(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			double lambda1 = ::Lerp(double(i + 1) / double(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
		}
		return r;
	}

	// magic
	void ToXYZ(double xyz[3]) const
	{
		xyz[0] = xyz[1] = xyz[2] = 0;
		for (int i = 0; i < nSpectralSamples; ++i) {
			xyz[0] += X.c[i] * c[i];
			xyz[1] += Y.c[i] * c[i];
			xyz[2] += Z.c[i] * c[i];
		}
		double scale = double(sampledLambdaEnd - sampledLambdaStart) /
			double(CIE_Y_integral * nSpectralSamples);
		xyz[0] *= scale;
		xyz[1] *= scale;
		xyz[2] *= scale;
	}

	// more magic
	double y() const
	{
		double yy = 0;
		for (int i = 0; i < nSpectralSamples; ++i)
			yy += Y.c[i] * c[i];
		return yy * double(sampledLambdaEnd - sampledLambdaStart) /
			double(CIE_Y_integral * nSpectralSamples);
	}

	// Compute XYZ matching functions then RGB to spectrum functions
	// with some weird shit and magic
	static void Init()
	{
		for (int i = 0; i < nSpectralSamples; ++i) {
			double wl0 = ::Lerp(double(i) / double(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			double wl1 = ::Lerp(double(i + 1) / double(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples,
				wl0, wl1);
			Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples,
				wl0, wl1);
			Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples,
				wl0, wl1);
		}
	}

private:
	// global X Y Z matching functions
	static SampledSpectrum X, Y, Z;

	// TODO: p325 RGB color
};


class RGBSpectrum : public CoefficientSpectrum<3>
{
public:
	RGBSpectrum(double v = 0) : CoefficientSpectrum<3>(v) { }
	RGBSpectrum(const CoefficientSpectrum<3>& v)
		: CoefficientSpectrum<3>(v) { }

	// converting RGB values to spectrum is trivial
	static RGBSpectrum FromRGB(const double rgb[3],
		SpectrumType type = SpectrumType::Reflectance) 
	{
		RGBSpectrum s;
		s.c[0] = rgb[0];
		s.c[1] = rgb[1];
		s.c[2] = rgb[2];
		return s;
	}

	// same as above
	void ToRGB(double* rgb) const 
	{
		rgb[0] = c[0];
		rgb[1] = c[1];
		rgb[2] = c[2];
	}

	// but why
	const RGBSpectrum& ToRGBSpectrum() const 
	
	{
		return *this;
	}

	void ToXYZ() 
	{
	// TODO: 
	}

	static RGBSpectrum FromXYZ(double *XYZ)
	{
		// TODO:
		return NULL;
	}

	void y()
	{
		// TODO:
	}

	// converts the spectrum to XYZand then to RGB.
	static RGBSpectrum FromSampled(const double* lambda, const double* v,
		int n) {

		if (!SpectrumSamplesSorted(lambda, v, n)) {
			// create a new copy to this data
			std::vector<double> slambda(&lambda[0], &lambda[n]);
			std::vector<double> sv(&v[0], &v[n]);
			// sort it
			SortSpectrumSamples(&slambda[0], &sv[0], n);
			return FromSampled(&slambda[0], &sv[0], n);
		}
		
		double xyz[3] = { 0, 0, 0 };
		for (int i = 0; i < nCIESamples; ++i) {
			double val = InterpolateSpectrumSamples(lambda, v, n,
				CIE_lambda[i]);
			xyz[0] += val * CIE_X[i];
			xyz[1] += val * CIE_Y[i];
			xyz[2] += val * CIE_Z[i];
		}
		double scale = double(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
			double(CIE_Y_integral * nCIESamples);
		xyz[0] *= scale;
		xyz[1] *= scale;
		xyz[2] *= scale;
		return FromXYZ(xyz);
	}


};