#pragma once

#include "core/texture.h"


// simply convert (u, v) to green and red
class UVTexture : public Texture<Spectrum> 
{
private:
	std::unique_ptr<TextureMapping2D> mapping;

public:
	Spectrum Evaluate(const SurfaceInteraction& si) const 
	{
		Vector2f dstdx, dstdy;
		Point2f st = mapping->Map(si, &dstdx, &dstdy);
		Float rgb[3] = { st[0] - std::floor(st[0]), st[1] - std::floor(st[1]), 0 };
		return Spectrum::FromRGB(rgb);
	}

};