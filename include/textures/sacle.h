#pragma once

#include "core/texture.h"

// takes two textures and returns the product of their values when evaluated
template <typename T1, typename T2>
class ScaleTexture : public Texture<T2> 
{
private:
	std::shared_ptr<Texture<T1>> tex1;
	std::shared_ptr<Texture<T2>> tex2;

public:
	ScaleTexture(const std::shared_ptr<Texture<T1>>& tex1,
		const std::shared_ptr<Texture<T2>>& tex2)
		: tex1(tex1), tex2(tex2) { }

	T2 Evaluate(const SurfaceInteraction& si) const 
	{
		return tex1->Evaluate(si) * tex2->Evaluate(si);
	}
};