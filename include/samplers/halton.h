#pragma once
#include "core/sampler.h"
#include "utlis/lowdiscrepancy.h"

static PBRT_CONSTEXPR int kMaxResolution = 128;

// HaltonSampler Utility Functions
static void extendedGCD(uint64_t a, uint64_t b, int64_t* x, int64_t* y);
static uint64_t multiplicativeInverse(int64_t a, int64_t n) 
{
    int64_t x, y;
    extendedGCD(a, n, &x, &y);
    return Mod(x, n);
}

static void extendedGCD(uint64_t a, uint64_t b, int64_t* x, int64_t* y)
{
    if (b == 0) {
        *x = 1;
        *y = 0;
        return;
    }
    int64_t d = a / b, xp, yp;
    extendedGCD(b, a % b, &xp, &yp);
    *x = yp;
    *y = xp - (d * yp);
}

// Sampler using Halton sequence to generate
// samples with good discrepancy 
class HaltonSampler : public GlobalSampler
{
private:
    // HaltonSampler Private Data
    static std::vector<uint16_t> radicalInversePermutations;
    Point2i baseScales, baseExponents;
    int sampleStride;
    int multInverse[2];
    mutable Point2i pixelForOffset = Point2i(std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max());
    mutable int64_t offsetForCurrentPixel;
    // Added after book publication: force all image samples to be at the
    // center of the pixel area.
    bool sampleAtPixelCenter;

public:
	HaltonSampler::HaltonSampler(int samplesPerPixel,
		const Bounds2i& sampleBounds);

    int64_t GetIndexForSample(int64_t sampleNum) const;

    Float SampleDimension(int64_t index, int dimension) const;

    std::unique_ptr<Sampler> Clone(int seed);
	
private:
	const uint16_t* PermutationForDimension(int dim) const 
	{
		return &radicalInversePermutations[PrimeSums[dim]];
	}


};
