#pragma once

#include "core/texture.h"

template <typename T> class DotsTexture : public Texture<T> 
{
private:
	std::unique_ptr<TextureMapping2D> mapping;
	std::shared_ptr<Texture<T>> outsideDot, insideDot;

public:
	DotsTexture(std::unique_ptr<TextureMapping2D> mapping,
		const std::shared_ptr<Texture<T>>& outsideDot,
		const std::shared_ptr<Texture<T>>& insideDot)
		: mapping(std::move(mapping)), outsideDot(outsideDot),
		insideDot(insideDot) { }

	T Evaluate(const SurfaceInteraction& si) const
	{
		Vector2f dstdx, dstdy;
		Point2f st = mapping->Map(si, &dstdx, &dstdy);
		int sCell = std::floor(st[0] + 0.5), tCell = std::floor(st[1] + 0.5);

		if (Noise(sCell + 0.5, tCell + 0.5) > 0) 
		{
			Float radius = 0.35;
			Float maxShift = 0.5 - radius;
			Float sCenter = sCell +
				maxShift * Noise(sCell + 1.5, tCell + 2.8);
			Float tCenter = tCell +
				maxShift * Noise(sCell + 4.5, tCell + 9.8);
			Vector2f dst = st - Point2f(sCenter, tCenter);
			if (dst.LengthSquared() < radius * radius)
				return insideDot->Evaluate(si);
		}
		return outsideDot->Evaluate(si);
	}
};