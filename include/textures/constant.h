#pragma once
#include "core/texture.h"

// a texture return the same value everywhere
template <typename T>
class ConstantTexture : public Texture<T>
{
private:
	T value;

public:
	ConstantTexture(const T& value) : value(value)
	{ }
	
	T Evaluate(const SurfaceInteraction&) const 
	{
		return value;
	}

};

