#pragma once
#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/interaction.h"
#include "core/transform.h"

class TextureMapping2D 
{
public:
	virtual Point2f Map(const SurfaceInteraction& si,
		Vector2f* dstdx, Vector2f* dstdy) const = 0;

	virtual ~TextureMapping2D() {}
};

class UVMapping2D : public TextureMapping2D 
{
private:
	const Float su, sv, du, dv;

public:
	UVMapping2D(Float su, Float sv, Float du, Float dv)
		: su(su), sv(sv), du(du), dv(dv) { }

	Point2f Map(const SurfaceInteraction& si,
		Vector2f* dstdx, Vector2f* dstdy) const;

};

class SphericalMapping2D : public TextureMapping2D 
{

private:
	const Transform WorldToTexture;

public:
	SphericalMapping2D(const Transform& WorldToTexture)
		: WorldToTexture(WorldToTexture) {}

	Point2f Map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const;
	
private:
	Point2f sphere(const Point3f& p) const;


};

class CylindricalMapping2D : public TextureMapping2D 
{
private:
	const Transform WorldToTexture;

public:
	CylindricalMapping2D(const Transform& WorldToTexture)
		: WorldToTexture(WorldToTexture) {}

	Point2f Map(const SurfaceInteraction& si,
		Vector2f* dstdx, Vector2f* dstdy) const;

private:
	Point2f cylinder(const Point3f& p) const;

};

class PlanarMapping2D : public TextureMapping2D 
{
private:
	const Vector3f vs, vt;
	const Float ds, dt;

public:
	PlanarMapping2D(const Vector3f& vs, const Vector3f& vt,
		Float ds = 0, Float dt = 0)
		: vs(vs), vt(vt), ds(ds), dt(dt) { }

	Point2f Map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const;

};

class TextureMapping3D 
{
public:
	virtual Point3f Map(const SurfaceInteraction& si,
		Vector3f* dpdx, Vector3f* dpdy) const = 0;

	virtual ~TextureMapping3D() {}
};

class TransformMapping3D : public TextureMapping3D 
{
private:
	const Transform WorldToTexture;

public:
	Point3f Map(const SurfaceInteraction& si,
		Vector3f* dpdx, Vector3f* dpdy) const;

};


template <typename T>
class Texture
{
public:
	virtual ~Texture() {}

	virtual T Evaluate(const SurfaceInteraction&) const = 0;

private:

};

// noises, perlin and fbm
static constexpr int NoisePermSize = 256;
static int NoisePerm[2 * NoisePermSize] = {
	151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140,
	36, 103, 30, 69, 142,
	// Remainder of the noise permutation table
	8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62,
	94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174,
	20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77,
	146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55,
	46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76,
	132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100,
	109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147,
	118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28,
	42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101,
	155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
	178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12,
	191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31,
	181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
	138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66,
	215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194,
	233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6,
	148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32,
	57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74,
	165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60,
	211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25,
	63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135,
	130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226,
	250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59,
	227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2,
	44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19,
	98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251,
	34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249,
	14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115,
	121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72,
	243, 141, 128, 195, 78, 66, 215, 61, 156, 180 };

inline Float NoiseWeight(Float t) {
	Float t3 = t * t * t;
	Float t4 = t3 * t;
	return 6 * t4 * t - 15 * t4 + 10 * t3;
}

inline Float Grad(int x, int y, int z, Float dx, Float dy, Float dz) 
{
	int h = NoisePerm[NoisePerm[NoisePerm[x] + y] + z];
	h &= 15;
	Float u = h < 8 || h == 12 || h == 13 ? dx : dy;
	Float v = h < 4 || h == 12 || h == 13 ? dy : dz;
	return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
}

inline Float SmoothStep(Float min, Float max, Float value) {
	Float v = Clamp((value - min) / (max - min), 0, 1);
	return v * v * (-2 * v + 3);
}

Float Noise(Float x, Float y = 0.5, Float z = 0.5);
Float Noise(const Point3f& p);

// TODO: something something music related 
Float FBm(const Point3f& p, const Vector3f& dpdx, const Vector3f& dpdy,
	Float omega, int maxOctaves);

Float Turbulence(const Point3f& p, const Vector3f& dpdx,
	const Vector3f& dpdy, Float omega, int maxOctaves);