#pragma once
#include "core/pbrt.h"
#include "core/camera.h"
#include "core/rng.h"


class Sampler
{
public:
	const int64_t samplesPerPixel;

protected:
	Point2i currentPixel;
	int64_t currentPixelSampleIndex;
	std::vector<int> samples1DArraySizes, samples2DArraySizes;
	std::vector<std::vector<Float>> sampleArray1D;
	std::vector<std::vector<Point2f>> sampleArray2D;

private:
	size_t array1DOffset, array2DOffset;

public:
	Sampler(int64_t samplesPerPixel);
	~Sampler() {}

	// start working on a pixel
	virtual void StartPixel(const Point2i& p);

	// next dimension of the current sample vector
	virtual Float Get1D() = 0;
	// next 2 dimensions of the current sample vector
	virtual Point2f Get2D() = 0;

	void Request1DArray(int n);
	void Request2DArray(int n);
	const Float* Get1DArray(int n);
	const Point2f* Get2DArray(int n);

	virtual int RoundCount(int n) const 
	{
		return n;
	}

	virtual bool StartNextSample();

	virtual std::unique_ptr<Sampler> Clone(int seed) = 0;

	// set the index of the sample in the current pixel
	// This method returns false once sampleNum is greater than or equal to samplesPerPixel
	virtual bool SetSampleNumber(int64_t sampleNum);

	CameraSample GetCameraSample(Point2i pRaster);
};

// sample based on each pixel
class PixelSampler : public Sampler 
{

protected:
	std::vector<std::vector<Float>> samples1D;
	std::vector<std::vector<Point2f>> samples2D;
	int current1DDimension = 0, current2DDimension = 0;
	RNG rng;

public:
	PixelSampler(int64_t samplesPerPixel, int nSampledDimensions);

	std::unique_ptr<Sampler> Clone(int seed);

	bool StartNextSample();

	bool SetSampleNumber(int64_t sampleNum);

	Float Get1D();


	
};

// sample based on the global scene then map to each pixel
class GlobalSampler : public Sampler
{
private:
	int dimension;
	int64_t intervalSampleIndex;
	static const int arrayStartDim = 5;
	int arrayEndDim;

public:
	GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) { }

	virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;

	virtual Float SampleDimension(int64_t index, int dimension) const = 0;

	void StartPixel(const Point2i& p);

	bool StartNextSample();
	bool SetSampleNumber(int64_t sampleNum);
	Float Get1D();
	Point2f Get2D();
};


// helper functions
void StratifiedSample1D(Float* samp, int nSamples, RNG& rng, bool jitter);

void StratifiedSample2D(Point2f* samp, int nx, int ny, RNG& rng, bool jitter);

template <typename T>
inline void Shuffle(T* samp, int count, int nDimensions, RNG& rng) 
{
	for (int i = 0; i < count; ++i) {
		int other = i + rng.UniformUInt32(count - i);
		for (int j = 0; j < nDimensions; ++j)
			std::swap(samp[nDimensions * i + j],
				samp[nDimensions * other + j]);
	}
}

void LatinHypercube(Float* samples, int nSamples, int nDim, RNG& rng);