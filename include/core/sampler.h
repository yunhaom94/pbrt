#pragma once
#include "core/pbrt.h"
#include "core/camera.h"

// TODO: p421, implement methods and move internal includes in .cpp file
class Sampler
{

public:
		Float samplesPerPixel;

public:
	Sampler();
	~Sampler();

	std::unique_ptr<Sampler> Clone(int seed)
	{
		return nullptr;
	}

	void StartPixel(Point2i pixel)
	{

	}

	bool StartNextSample()
	{

	}

	CameraSample GetCameraSample(Point2i pixel)
	{
		return *(new CameraSample);
	}

private:

};

Sampler::Sampler()
{
}

Sampler::~Sampler()
{
}