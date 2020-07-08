#include "core/spectrum.h"

SampledSpectrum SampledSpectrum::X;
SampledSpectrum SampledSpectrum::Y;
SampledSpectrum SampledSpectrum::Z;
SampledSpectrum SampledSpectrum::rgbRefl2SpectWhite;
SampledSpectrum SampledSpectrum::rgbRefl2SpectCyan;
SampledSpectrum SampledSpectrum::rgbRefl2SpectMagenta;
SampledSpectrum SampledSpectrum::rgbRefl2SpectYellow;
SampledSpectrum SampledSpectrum::rgbRefl2SpectRed;
SampledSpectrum SampledSpectrum::rgbRefl2SpectGreen;
SampledSpectrum SampledSpectrum::rgbRefl2SpectBlue;
SampledSpectrum SampledSpectrum::rgbIllum2SpectWhite;
SampledSpectrum SampledSpectrum::rgbIllum2SpectCyan;
SampledSpectrum SampledSpectrum::rgbIllum2SpectMagenta;
SampledSpectrum SampledSpectrum::rgbIllum2SpectYellow;
SampledSpectrum SampledSpectrum::rgbIllum2SpectRed;
SampledSpectrum SampledSpectrum::rgbIllum2SpectGreen;
SampledSpectrum SampledSpectrum::rgbIllum2SpectBlue;


void Blackbody(const Float* lambda, int n, Float T, Float* Le)
{
    const Float c = 299792458;
    const Float h = 6.62606957e-34;
    const Float kb = 1.3806488e-23;
    for (int i = 0; i < n; ++i) 
    {
        Float l = lambda[i] * 1e-9;
        Float lambda5 = (l * l) * (l * l) * l;
        Le[i] = (2 * h * c * c) /
            (lambda5 * (std::exp((h * c) / (l * kb * T)) - 1));
    }
}

void BlackbodyNormalized(const Float* lambda, int n, Float T, Float* Le)
{
    Blackbody(lambda, n, T, Le);
    Float lambdaMax = 2.8977721e-3 / T * 1e9;
    Float maxL;
    Blackbody(&lambdaMax, 1, T, &maxL);
    for (int i = 0; i < n; ++i)
        Le[i] /= maxL;
}

bool SpectrumSamplesSorted(const Float* lambda, const Float* v, int n)
{
	Float last_l = lambda[0];

	for (int i = 0; i < n; i++)
	{
		Float l = lambda[i];
		// check if wavelength is not sorted
		if (l < last_l)
			return false;

		last_l = l;
	}

	return true;
}

void SortSpectrumSamples(Float* lambda, Float* v, int n)
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

// For each such segment, it computes the average
// value over its range, scales the average by the wavelength range the segment covers, and
// accumulates a sum of these values.
Float AverageSpectrumSamples(const Float* lambda, const Float* vals, int n, Float lambdaStart, Float lambdaEnd)
{
	//boundary check
	if (lambdaEnd <= lambda[0]) 
		return vals[0];
	if (lambdaStart >= lambda[n - 1])
		return vals[n - 1];
	if (n == 1)
		return vals[0];

	Float sum = 0;

	// Add contributions of constant segments before/after samples
	if (lambdaStart < lambda[0])
		sum += vals[0] * (lambda[0] - lambdaStart);
	if (lambdaEnd > lambda[n - 1])
		sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);

	int i = 0;
	while (lambdaStart > lambda[i + 1])
		++i;

	auto interp = [lambda, vals](Float w, int i)
	{
		return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]),
			vals[i], vals[i + 1]);
	};

	for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) 
	{
		Float segLambdaStart = std::max(lambdaStart, lambda[i]);
		Float segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
		sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
			(segLambdaEnd - segLambdaStart);
	}

	return sum / (lambdaEnd - lambdaStart);
}

Float InterpolateSpectrumSamples(const Float* lambda, const Float* vals,
    int n, Float l) 
{
    if (l <= lambda[0]) return vals[0];
    if (l >= lambda[n - 1]) return vals[n - 1];
    int offset = FindInterval(n,
        [&](int index) { return lambda[index] <= l; });
    Float t = (l - lambda[offset]) / (lambda[offset + 1] - lambda[offset]);
    return Lerp(t, vals[offset], vals[offset + 1]);
}

SampledSpectrum SampledSpectrum::FromRGB(const Float rgb[3], SpectrumType type)
{
    SampledSpectrum r;

    if (type == SpectrumType::Reflectance) {
        // Convert reflectance spectrum to RGB
        if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]) {
            // Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum
            r += rgb[0] * rgbRefl2SpectWhite;
                
            if (rgb[1] <= rgb[2]) {
                r += (rgb[1] - rgb[0]) * rgbRefl2SpectCyan;
                r += (rgb[2] - rgb[1]) * rgbRefl2SpectBlue;
            }
            else {
                r += (rgb[2] - rgb[0]) * rgbRefl2SpectCyan;
                r += (rgb[1] - rgb[2]) * rgbRefl2SpectGreen;
            }
        }
        else if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]) {
            // Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
            r += rgb[1] * rgbRefl2SpectWhite;
            if (rgb[0] <= rgb[2]) {
                r += (rgb[0] - rgb[1]) * rgbRefl2SpectMagenta;
                r += (rgb[2] - rgb[0]) * rgbRefl2SpectBlue;
            }
            else {
                r += (rgb[2] - rgb[1]) * rgbRefl2SpectMagenta;
                r += (rgb[0] - rgb[2]) * rgbRefl2SpectRed;
            }
        }
        else {
            // Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
            r += rgb[2] * rgbRefl2SpectWhite;
            if (rgb[0] <= rgb[1]) {
                r += (rgb[0] - rgb[2]) * rgbRefl2SpectYellow;
                r += (rgb[1] - rgb[0]) * rgbRefl2SpectGreen;
            }
            else {
                r += (rgb[1] - rgb[2]) * rgbRefl2SpectYellow;
                r += (rgb[0] - rgb[1]) * rgbRefl2SpectRed;
            }
        }
        r *= .94;
    }
    else {
        // Convert illuminant spectrum to RGB
        if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]) {
            // Compute illuminant _SampledSpectrum_ with _rgb[0]_ as minimum
            r += rgb[0] * rgbIllum2SpectWhite;
            if (rgb[1] <= rgb[2]) {
                r += (rgb[1] - rgb[0]) * rgbIllum2SpectCyan;
                r += (rgb[2] - rgb[1]) * rgbIllum2SpectBlue;
            }
            else {
                r += (rgb[2] - rgb[0]) * rgbIllum2SpectCyan;
                r += (rgb[1] - rgb[2]) * rgbIllum2SpectGreen;
            }
        }
        else if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]) {
            // Compute illuminant _SampledSpectrum_ with _rgb[1]_ as minimum
            r += rgb[1] * rgbIllum2SpectWhite;
            if (rgb[0] <= rgb[2]) {
                r += (rgb[0] - rgb[1]) * rgbIllum2SpectMagenta;
                r += (rgb[2] - rgb[0]) * rgbIllum2SpectBlue;
            }
            else {
                r += (rgb[2] - rgb[1]) * rgbIllum2SpectMagenta;
                r += (rgb[0] - rgb[2]) * rgbIllum2SpectRed;
            }
        }
        else {
            // Compute illuminant _SampledSpectrum_ with _rgb[2]_ as minimum
            r += rgb[2] * rgbIllum2SpectWhite;
            if (rgb[0] <= rgb[1]) {
                r += (rgb[0] - rgb[2]) * rgbIllum2SpectYellow;
                r += (rgb[1] - rgb[0]) * rgbIllum2SpectGreen;
            }
            else {
                r += (rgb[1] - rgb[2]) * rgbIllum2SpectYellow;
                r += (rgb[0] - rgb[1]) * rgbIllum2SpectRed;
            }
        }
        r *= 0.86445;
    }

    return r.Clamp();

}

SampledSpectrum::SampledSpectrum(const RGBSpectrum& r, SpectrumType t)
{
    Float rgb[3];
    r.ToRGB(rgb);
    *this = SampledSpectrum::FromRGB(rgb, t);
}

RGBSpectrum SampledSpectrum::ToRGBSpectrum() const
{


    return RGBSpectrum();
}
