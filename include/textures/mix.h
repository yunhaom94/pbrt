#pragma once

#include "core/texture.h"

// lerp between 2 texture
template <typename T>
class MixTexture : public Texture<T> 
{
private:
	std::shared_ptr<Texture<T>> tex1, tex2;
	std::shared_ptr<Texture<Float>> amount;

public:
	MixTexture(const std::shared_ptr<Texture<T>>& tex1,
		const std::shared_ptr<Texture<T>>& tex2,
		const std::shared_ptr<Texture<Float>>& amount)
		: tex1(tex1), tex2(tex2), amount(amount) { }

	T Evaluate(const SurfaceInteraction& si) const
	{
		T t1 = tex1->Evaluate(si), t2 = tex2->Evaluate(si);
		Float amt = amount->Evaluate(si);
		return (1 - amt) * t1 + amt * t2;
	}

};