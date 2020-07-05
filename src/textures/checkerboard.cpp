#include "textures/checkerboard.h"

template<typename T>
inline T Checkerboard2DTexture<T>::Evaluate(const SurfaceInteraction& si) const
{
	Vector2f dstdx, dstdy;
	Point2f st = mapping->Map(si, &dstdx, &dstdy);
	if (aaMethod == AAMethod::None) {
		if (((int)std::floor(st[0]) + (int)std::floor(st[1])) % 2 == 0)
			return tex1->Evaluate(si);
		return tex2->Evaluate(si);
	}
	else // this is only a bx filter
	{
		Float ds = std::max(std::abs(dstdx[0]), std::abs(dstdy[0]));
		Float dt = std::max(std::abs(dstdx[1]), std::abs(dstdy[1]));
		Float s0 = st[0] - ds, s1 = st[0] + ds;
		Float t0 = st[1] - dt, t1 = st[1] + dt;
		if (std::floor(s0) == std::floor(s1) &&
			std::floor(t0) == std::floor(t1)) {
			if (((int)std::floor(st[0]) + (int)std::floor(st[1])) % 2 == 0)
				return tex1->Evaluate(si);
			return tex2->Evaluate(si);
		}
		auto bumpInt = [](Float x) {
			return (int)std::floor(x / 2) +
				2 * std::max(x / 2 - (int)std::floor(x / 2) - (Float)0.5,
					(Float)0); };
		Float sint = (bumpInt(s1) - bumpInt(s0)) / (2 * ds);
		Float tint = (bumpInt(t1) - bumpInt(t0)) / (2 * dt);
		Float area2 = sint + tint - 2 * sint * tint;
		if (ds > 1 || dt > 1)
			area2 = .5f;
		return (1 - area2) * tex1->Evaluate(si) +
			area2 * tex2->Evaluate(si);
	}
}
