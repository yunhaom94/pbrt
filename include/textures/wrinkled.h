#pragma once

#include "core/texture.h"

template <typename T>
class WrinkledTexture : public Texture<T>
{
private:
	std::unique_ptr<TextureMapping3D> mapping;
	const Float omega;
	const int octaves;

public:
	WrinkledTexture(std::unique_ptr<TextureMapping3D> mapping, int octaves,
		Float omega)
		: mapping(std::move(mapping)), omega(omega), octaves(octaves) { }

	T Evaluate(const SurfaceInteraction& si) const
	{
		Vector3f dpdx, dpdy;
		Point3f P = mapping->Map(si, &dpdx, &dpdy);
		return Turbulence(P, dpdx, dpdy, omega, octaves);
	}

};