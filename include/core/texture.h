#pragma once
#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/interaction.h"

// TODO: p614
template <typename T>
class Texture
{
public:
	Texture() {}
	~Texture() {}

	T Evaluate(SurfaceInteraction s) { return T(0); }

private:

};

