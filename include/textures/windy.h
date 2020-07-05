#pragma once

#include "core/texture.h"

// windy water surface
template <typename T> class WindyTexture : public Texture<T> 
{
private:
	std::unique_ptr<TextureMapping3D> mapping;

public:
	WindyTexture(std::unique_ptr<TextureMapping3D> mapping)
		: mapping(std::move(mapping)) { }

	Evaluate(const SurfaceInteraction& si) const {
		Vector3f dpdx, dpdy;
		Point3f P = mapping->Map(si, &dpdx, &dpdy);
		Float windStrength = FBm(.1f * P, .1f * dpdx, .1f * dpdy, .5, 3);
		Float waveHeight = FBm(P, dpdx, dpdy, .5, 6);
		return std::abs(windStrength) * waveHeight;
	}
};