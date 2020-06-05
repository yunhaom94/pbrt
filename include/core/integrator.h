#pragma once
#include "core/pbrt.h"

// This is the class with all the BSDF algorithms
// it is able to actually render a image
class Integrator
{
public:
	virtual void Render(const Scene& scene) = 0;

private:

};

