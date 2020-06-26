#pragma once
#include "core/pbrt.h"
#include "core/sampler.h"



// Split each pixel into different non-overlapping blocks
// Create one sample for each block
class StratifiedSampler : public PixelSampler 
{
private:
	const int xPixelSamples, yPixelSamples;
	const bool jitterSamples;

public:
	StratifiedSampler(int xPixelSamples, int yPixelSamples,
		bool jitterSamples, int nSampledDimensions);
	void StartPixel(const Point2i& p);


}
};