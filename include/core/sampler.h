#pragma once
#include "core/pbrt.h"
#include "core/camera.h"
#include "core/sampling.h"
#include "utlis/rng.h"



class Sampler
{
	/* How sampler stores data?
	* In a n-dimension width 2D array
	* Dimensions/Samples:       s1,  s2,  s3, ...
	* x (x coord in pixel) :   [0.2, 0,4 ,0.5 ...] 
	* y (y coord in pixel) :   [0.3, 0.2, 0.1 ...]
	* t (time) :			   [0.8, 0.4, 0.6 ...]
	* u :					   [0.3, 0.9, 0.7 ...]
	* v :					   [0.1, 0.5, 0.3 ...]
	* More dimensions...	   [...]
	*/

public:
	const int64_t samplesPerPixel;

protected:
	Point2i currentPixel;
	int64_t currentPixelSampleIndex;

	// sample values for each dimension stored in arrays
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

	// compute next sample, return false if more
	// sample then samples per pixel
	virtual bool StartNextSample();

	// returns the sample value for next 
	// dimension of the current sample vector.
	// kinda of like an iterator,
	// For example, current at sample value for x, call it again return y
	virtual Float Get1D() = 0;

	// returns the sample values for the next two dimensions
	virtual Point2f Get2D() = 0;

	// If the entire array of samples of some dimension are needed, 
	// they must be requested before rendering begins.
	void Request1DArray(int n);
	void Request2DArray(int n);

	// return pointer to the start of a previously 
	// requested array of samples for the current dimension.
	const Float* Get1DArray(int n);
	const Point2f* Get2DArray(int n);

	// find the optimal sample count for this sampler
	// given the desired sample count
	virtual int RoundCount(int n) const 
	{
		return n;
	}

	virtual std::unique_ptr<Sampler> Clone(int seed) = 0;

	// set the index of the sample in the current pixel
	// This method returns false once sampleNum is greater than or equal to samplesPerPixel
	virtual bool SetSampleNumber(int64_t sampleNum);

	// initializes a CameraSample for a given pixel.
	CameraSample GetCameraSample(Point2i pRaster);
};

// give sample points based on each pixel
class PixelSampler : public Sampler 
{

protected:
	// sample1D[dim][pixelSample]
	std::vector<std::vector<Float>> samples1D;
	std::vector<std::vector<Point2f>> samples2D;
	int current1DDimension = 0, current2DDimension = 0;
	RNG rng;

public:
	PixelSampler(int64_t samplesPerPixel, int nSampledDimensions);

	bool StartNextSample();

	bool SetSampleNumber(int64_t sampleNum);

	Float Get1D();
	Point2f Get2D();
};

// sample on the entire imagine plane then map each sample point to each pixel
class GlobalSampler : public Sampler
{
private:
	int dimension;
	int64_t intervalSampleIndex;
	static const int arrayStartDim = 5;
	int arrayEndDim;

public:
	GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) { }

	// given index of sample point for current pixel, 
	// get the index for the global sample vector
	virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;

	// returns the sample value for the given dimension of
	// the indexth sample vector in the sequence.
	virtual Float SampleDimension(int64_t index, int dimension) const = 0;

	void StartPixel(const Point2i& p);

	bool StartNextSample();
	bool SetSampleNumber(int64_t sampleNum);
	Float Get1D();
	Point2f Get2D();
};


