#pragma once

#include "core/texture.h"

enum class AAMethod { None, ClosedForm };

template <typename T> 
class Checkerboard2DTexture : public Texture<T>
{
private:
	std::unique_ptr<TextureMapping2D> mapping;
	const std::shared_ptr<Texture<T>> tex1, tex2;
	const AAMethod aaMethod;

public:
	Checkerboard2DTexture(std::unique_ptr<TextureMapping2D> mapping,
		const std::shared_ptr<Texture<T>>& tex1,
		const std::shared_ptr<Texture<T>>& tex2, AAMethod aaMethod)
		: mapping(std::move(mapping)), tex1(tex1), tex2(tex2),
		aaMethod(aaMethod) { }

	T Evaluate(const SurfaceInteraction& si) const;

};

template <typename T>
class Checkerboard3DTexture : public Texture<T> 
{
private:
	std::unique_ptr<TextureMapping3D> mapping;
	std::shared_ptr<Texture<T>> tex1, tex2;

public:
	Checkerboard3DTexture(std::unique_ptr<TextureMapping3D> mapping,
		const std::shared_ptr<Texture<T>>& tex1,
		const std::shared_ptr<Texture<T>>& tex2)
		: mapping(std::move(mapping)), tex1(tex1), tex2(tex2) { }

	T Evaluate(const SurfaceInteraction& si) const {
		Vector3f dpdx, dpdy;
		Point3f p = mapping->Map(si, &dpdx, &dpdy);
		if (((int)std::floor(p.x) + (int)std::floor(p.y) +
			(int)std::floor(p.z)) % 2 == 0)
			return tex1->Evaluate(si);
		else
			return tex2->Evaluate(si);
	}

};