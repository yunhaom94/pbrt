
#include "core/transform.h"
#include "core/interaction.h"
#include "core/ray.h"
#include "core/bounding_boxes.h"



Transform Inverse(const Transform& t)
{
	return Transform(t.mInv, t.m);
}

Transform Transpose(const Transform& t)
{
	return Transform(t.m.transpose(), t.mInv.transpose());
}

Transform Translate(const Vector3f& delta)
{
	Matrix4x4 m;
	m << 1, 0, 0, delta.x(),
		 0, 1, 0, delta.y(),
		 0, 0, 1, delta.z(),
		 0, 0, 0, 1;

	Matrix4x4 minv;
	minv << 1, 0, 0, -delta.x(),
		    0, 1, 0, -delta.y(),
		    0, 0, 1, -delta.z(),
		    0, 0, 0, 1;

	return Transform(m, minv);
}

Transform Scale(Float x, Float y, Float z)
{
	Matrix4x4 m;
	m << x, 0, 0, 0,
		 0, y, 0, 0,
		 0, 0, z, 0,
		 0, 0, 0, 1;

	Matrix4x4 minv;
	m << 1.0 / x, 0, 0, 0,
		 0, 1.0 / y, 0, 0,
		 0, 0, 1.0 / z, 0,
		 0, 0, 0, 1;

	return Transform(m, minv);
}

Transform RotateX(Float theta)
{
	Matrix4x4 m;
	m << 1, 0,               0,                0,
		 0, std::cos(theta), -std::sin(theta), 0,
		 0, std::sin(theta), std::cos(theta),  0,
		 0, 0,               0,                1;

	return Transform(m, m.transpose());
}

Transform RotateY(Float theta)
{
	Matrix4x4 m;
	m << std::cos(theta),  0, std::sin(theta), 0,
		 0,                1, 0,               0,
		 -std::sin(theta), 0, std::cos(theta), 0,
		 0,                0, 0,               1;


	return Transform(m, m.transpose());
}

Transform RotateZ(Float theta)
{
	Matrix4x4 m;
	m << std::cos(theta), -std::sin(theta), 0, 0,
		 std::sin(theta), std::cos(theta),  0, 0,
		 0,               0,                1, 0,
		 0,               0,                0, 1;


	return Transform(m, m.transpose());
}

Transform Rotate(Float theta, const Vector3f& axis)
{

	Vector3f a = axis.normalized();
	Float sinTheta = std::sin(theta);
	Float cosTheta = std::cos(theta);
	Matrix4x4 m = Matrix4x4::Zero();

	m.row(0)[0] = a.x() * a.x() + (1 - a.x() * a.x()) * cosTheta;
	m.row(0)[1] = a.x() * a.y() * (1 - cosTheta) - a.z() * sinTheta;
	m.row(0)[2] = a.x() * a.z() * (1 - cosTheta) + a.y() * sinTheta;
	m.row(0)[3] = 0;

	m.row(1)[0] = a.x() * a.y() * (1 - cosTheta) + a.z() * sinTheta;
	m.row(1)[1] = a.y() * a.y() + (1 - a.y() * a.y()) * cosTheta;
	m.row(1)[2] = a.y() * a.z() * (1 - cosTheta) - a.x() * sinTheta;
	m.row(1)[3] = 0;

	m.row(1)[0] = a.x() * a.z() * (1 - cosTheta) - a.y() * sinTheta;
	m.row(1)[1] = a.y() * a.z() * (1 - cosTheta) + a.x() * sinTheta;
	m.row(1)[2] = a.z() * a.z() + (1 - a.z() * a.z()) * cosTheta;
	m.row(1)[3] = 0;

	return Transform();
}

Transform LookAt(const Vector3f& pos, const Vector3f& look, const Vector3f& up)
{

	Vector3f dir = (look - pos).normalized();
	Vector3f left = (up.normalized().cross(dir)).normalized();
	Vector3f newUp = dir.cross(left);

	Matrix4x4 cameraToWorld;

	cameraToWorld << 
		left.x(), newUp.x(), dir.x(), pos.x(),
		left.y(), newUp.y(), dir.y(), pos.y(),
		left.z(), newUp.z(), dir.z(), pos.z(),
		0,        0,         0,       1;

	return Transform(cameraToWorld.inverse(), cameraToWorld);
}


Ray Transform::operator()(const Ray& r) const
{
	Vector3f oError;
	Point3f o3 = (*this)(r.o);
	Vector3f d3 = (*this)(r.d);

	// TODO: Offset ray origin to edge of error bounds and compute tMax p233

	return Ray(o3, d3, r.tMax, r.time, r.medium);
}

Ray Transform::operator()(const Ray& r, Vector3f* oError, Vector3f* dError) const
{
	Point3f o = (*this)(r.o, oError);
	Vector3f d = (*this)(r.d, dError);

	// TODO: Offset ray origin to edge of error bounds and compute tMax p233

	return Ray(o, d, r.tMax, r.time, r.medium);
}

Ray Transform::operator()(const Ray& r, const Vector3f& oErrorIn, const Vector3f& dErrorIn, Vector3f* oErrorOut, Vector3f* dErrorOut) const
{
	Point3f o = (*this)(r.o, oErrorIn, oErrorOut);
	Vector3f d = (*this)(r.d, dErrorIn, dErrorOut);

	// TODO: Offset ray origin to edge of error bounds and compute tMax p233

	return Ray(o, d, r.tMax, r.time, r.medium);
}


Bounds3f Transform::operator()(const Bounds3f& b) const
{
	const Transform& M = *this;
	Bounds3f ret(M(Point3f(b.pMin.x(), b.pMin.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMin.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMin.x(), b.pMax.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMin.x(), b.pMin.y(), b.pMax.z())));
	ret = Union(ret, M(Point3f(b.pMin.x(), b.pMax.y(), b.pMax.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMax.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMin.y(), b.pMax.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMax.y(), b.pMax.z())));
	return ret;


}

Transform Transform::operator*(const Transform& t2) const {
	return Transform(m * t2.m, mInv * t2.mInv);
}

bool Transform::HasScale() const
{
	Float la2 = (*this)(Vector3f(1, 0, 0)).squaredNorm();
	Float lb2 = (*this)(Vector3f(0, 1, 0)).squaredNorm();
	Float lc2 = (*this)(Vector3f(0, 0, 1)).squaredNorm();
#define NOT_ONE(x) ((x) < 0.999 || (x) > 1.001)
	return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
}

bool Transform::SwapsHandedness() const
{
	Eigen::Matrix3d mtemp;

	mtemp.row(0) = m.row(0).head<3>();
	mtemp.row(1) = m.row(1).head<3>();
	mtemp.row(2) = m.row(2).head<3>();

	return mtemp.determinant() < 0;
}

SurfaceInteraction Transform::operator()(const SurfaceInteraction& si) const
{
	SurfaceInteraction ret;
	// TODO: Transform p and pError in SurfaceInteraction p229

	// TODO: maybe write a clone function
	ret.wo = (*this)(si.wo);
	ret.n = (*this)(si.n);
	ret.time = si.time;
	ret.bsdf = si.bsdf;
	ret.surfaceInteraction = si.surfaceInteraction;
	ret.mediumInterface = si.mediumInterface;
	ret.uv = si.uv;
	ret.dpdu = (*this)(si.dpdu);
	ret.dpdv = (*this)(si.dpdv);
	ret.dndu = (*this)(si.dndu);
	ret.dndv = (*this)(si.dndv);
	ret.shape = si.shape;
	ret.shading.n = (*this)(si.shading.n).normalized();
	ret.shading.dpdu = (*this)(si.shading.dpdu);
	ret.shading.dpdv = (*this)(si.shading.dpdv);
	ret.shading.dndu = (*this)(si.shading.dndu);
	ret.shading.dndv = (*this)(si.shading.dndv);

	return ret;
}
