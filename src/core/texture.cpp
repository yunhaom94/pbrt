#include "core/texture.h"
#include "utlis/geometry.h"

Point2f UVMapping2D::Map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const
{
	*dstdx = Vector2f(su * si.dudx, sv * si.dvdx);
	*dstdy = Vector2f(su * si.dudy, sv * si.dvdy);

	return Point2f(su * si.uv[0] + du,
		sv * si.uv[1] + dv);
}

Point2f SphericalMapping2D::Map(const SurfaceInteraction& si,
	Vector2f* dstdx, Vector2f* dstdy) const {
	Point2f st = sphere(si.p);
	const Float delta = 0.1;

	Point2f stDeltaX = sphere(si.p + delta * si.dpdx);
	*dstdx = (stDeltaX - st) / delta;
	Point2f stDeltaY = sphere(si.p + delta * si.dpdy);
	*dstdy = (stDeltaY - st) / delta;

	if ((*dstdx)[1] > .5) 
		(*dstdx)[1] = 1 - (*dstdx)[1];
	else if
		((*dstdx)[1] < -.5f) (*dstdx)[1] = -((*dstdx)[1] + 1);
	
	if ((*dstdy)[1] > .5) 
		(*dstdy)[1] = 1 - (*dstdy)[1];
	else if 
		((*dstdy)[1] < -.5f) (*dstdy)[1] = -((*dstdy)[1] + 1);

	return st;
}

Point2f SphericalMapping2D::sphere(const Point3f& p) const
{
	Vector3f vec = (WorldToTexture(p) - Point3f(0, 0, 0)).normalized();
	Float theta = SphericalTheta(vec), phi = SphericalPhi(vec);
	return Point2f(theta * InvPi, phi * Inv2Pi);
}

Point2f CylindricalMapping2D::Map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const
{
	Point2f st = cylinder(si.p);
	// Compute texture coordinate differentials for cylinder $(u,v)$ mapping
	const Float delta = 0.01;
	Point2f stDeltaX = cylinder(si.p + delta * si.dpdx);
	*dstdx = (stDeltaX - st) / delta;
	if ((*dstdx)[1] > 0.5)
		(*dstdx)[1] = 1.0 - (*dstdx)[1];
	else if ((*dstdx)[1] < -0.5)
		(*dstdx)[1] = -((*dstdx)[1] + 1);
	Point2f stDeltaY = cylinder(si.p + delta * si.dpdy);
	*dstdy = (stDeltaY - st) / delta;
	if ((*dstdy)[1] > .5)
		(*dstdy)[1] = 1.0 - (*dstdy)[1];
	else if ((*dstdy)[1] < -0.5)
		(*dstdy)[1] = -((*dstdy)[1] + 1);
	return st;
}

Point2f CylindricalMapping2D::cylinder(const Point3f& p) const
{
	Vector3f vec = (WorldToTexture(p) - Point3f(0, 0, 0)).normalized();
	return Point2f((Pi + std::atan2(vec.y(), vec.x())) * Inv2Pi,
		vec.z());
}

Point2f PlanarMapping2D::Map(const SurfaceInteraction& si,
	Vector2f* dstdx, Vector2f* dstdy) const 
{
	Vector3f vec(si.p);
	*dstdx = Vector2f(si.dpdx.dot(vs), si.dpdx.dot(vt));
	*dstdy = Vector2f(si.dpdy.dot(vs), si.dpdy.dot(vt));
	return Point2f(ds + vec.dot(vs), dt + vec.dot(vt));
}

Point3f TransformMapping3D::Map(const SurfaceInteraction& si, Vector3f* dpdx, Vector3f* dpdy) const
{
	*dpdx = WorldToTexture(si.dpdx);
	*dpdy = WorldToTexture(si.dpdy);
	return WorldToTexture(si.p);
}

Float Noise(Float x, Float y, Float z)
{

	int ix = std::floor(x), iy = std::floor(y), iz = std::floor(z);
	Float dx = x - ix, dy = y - iy, dz = z - iz;

	ix &= NoisePermSize - 1;
	iy &= NoisePermSize - 1;
	iz &= NoisePermSize - 1;
	Float w000 = Grad(ix, iy, iz, dx, dy, dz);
	Float w100 = Grad(ix + 1, iy, iz, dx - 1, dy, dz);
	Float w010 = Grad(ix, iy + 1, iz, dx, dy - 1, dz);
	Float w110 = Grad(ix + 1, iy + 1, iz, dx - 1, dy - 1, dz);
	Float w001 = Grad(ix, iy, iz + 1, dx, dy, dz - 1);
	Float w101 = Grad(ix + 1, iy, iz + 1, dx - 1, dy, dz - 1);
	Float w011 = Grad(ix, iy + 1, iz + 1, dx, dy - 1, dz - 1);
	Float w111 = Grad(ix + 1, iy + 1, iz + 1, dx - 1, dy - 1, dz - 1);

	Float wx = NoiseWeight(dx), wy = NoiseWeight(dy), wz = NoiseWeight(dz);
	Float x00 = Lerp(wx, w000, w100);
	Float x10 = Lerp(wx, w010, w110);
	Float x01 = Lerp(wx, w001, w101);
	Float x11 = Lerp(wx, w011, w111);
	Float y0 = Lerp(wy, x00, x10);
	Float y1 = Lerp(wy, x01, x11);
	return Lerp(wz, y0, y1);

	return Float();
}

Float Noise(const Point3f& p)
{
	return Noise(p.x(), p.y(), p.z());
}

Float FBm(const Point3f& p, const Vector3f& dpdx, const Vector3f& dpdy, Float omega, int maxOctaves)
{
	Float len2 = std::max(dpdx.squaredNorm(), dpdy.squaredNorm());
	Float n = Clamp(-1 - .5f * Log2(len2), 0, maxOctaves);
	int nInt = std::floor(n);
	Float sum = 0, lambda = 1, o = 1;

	for (int i = 0; i < nInt; ++i) {
		sum += o * Noise(lambda * p);
		lambda *= 1.99f;
		o *= omega;
	}
	Float nPartial = n - nInt;
	sum += o * SmoothStep(.3f, .7f, nPartial) * Noise(lambda * p);
		return sum;
}

Float Turbulence(const Point3f& p, const Vector3f& dpdx, const Vector3f& dpdy, Float omega, int maxOctaves)
{
	Float len2 = std::max(dpdx.squaredNorm(), dpdy.squaredNorm());
	Float n = Clamp(-1 - .5f * Log2(len2), 0, maxOctaves);
	int nInt = std::floor(n);

	Float sum = 0, lambda = 1, o = 1;
	for (int i = 0; i < nInt; ++i) {
		sum += o * std::abs(Noise(lambda * p));
		lambda *= 1.99f;
		o *= omega;
	}

	Float nPartial = n - nInt;
	sum += o * Lerp(SmoothStep(.3f, .7f, nPartial),
		0.2, std::abs(Noise(lambda * p)));
	for (int i = nInt; i < maxOctaves; ++i) {
		sum += o * 0.2f;
		o *= omega;
	}
	return Float();
}

