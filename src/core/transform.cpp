
#include "core/transform.h"
#include "core/interaction.h"
#include "core/ray.h"
#include "core/bounding_boxes.h"



Transform Transpose(const Transform& t)
{
	return Transform(t.m.transpose(), t.mInv.transpose());
}

Transform Translate(const Eigen::Vector3d& delta)
{
	Eigen::Matrix4d m;
	m << 1, 0, 0, delta.x(),
		0, 1, 0, delta.y(),
		0, 0, 1, delta.z(),
		0, 0, 0, 1;

	Eigen::Matrix4d minv;
	minv << 1, 0, 0, -delta.x(),
		0, 1, 0, -delta.y(),
		0, 0, 1, -delta.z(),
		0, 0, 0, 1;



	return Transform(m, minv);
}

Transform Scale(double x, double y, double z)
{
	Eigen::Matrix4d m;
	m << x, 0, 0, 0,
		0, y, 0, 0,
		0, 0, z, 0,
		0, 0, 0, 1;

	Eigen::Matrix4d minv;
	m << 1/x, 0, 0, 0,
		0, 1/y, 0, 0,
		0, 0, 1/z, 0,
		0, 0, 0, 1;

	return Transform(m, minv);
}

Transform RotateX(double theta)
{
	Eigen::Matrix4d m;
	m << 1, 0, 0, 0,
		0, cos(theta), -sin(theta), 0,
		0, sin(theta), cos(theta), 0,
		0, 0, 0, 1;


	return Transform(m, m.transpose());
}

Transform RotateY(double theta)
{
	Eigen::Matrix4d m;
	m << cos(theta), 0, sin(theta), 0,
		0, 1, 0, 0,
		-sin(theta), 0, cos(theta), 0,
		0, 0, 0, 1;


	return Transform(m, m.transpose());
}

Transform RotateZ(double theta)
{
	Eigen::Matrix4d m;
	m << cos(theta), -sin(theta), 0, 0,
		sin(theta), cos(theta), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;


	return Transform(m, m.transpose());
}

Transform LookAt(const Eigen::Vector3d& pos, const Eigen::Vector3d& look, const Eigen::Vector3d& up)
{

	Eigen::Vector3d dir = (look - pos).normalized();
	Eigen::Vector3d left = (up.cross(dir)).normalized();
	Eigen::Vector3d newUp = dir.cross(left);

	Eigen::Matrix4d cameraToWorld;

	cameraToWorld << 
		left.x(), newUp.x(), dir.x(), pos.x(),
		left.y(), newUp.y(), dir.y(), pos.y(),
		left.z(), newUp.z(), dir.z(), pos.z(),
		0, 0, 0, 1;

	return Transform(cameraToWorld.inverse(), cameraToWorld);
}

Eigen::Vector3d Transform::operator()(const Eigen::Vector3d v) const
{
	Eigen::Vector3d ret = (Eigen::Vector4d(v.x(), v.y(), v.z(), 0).transpose() * m).head<3>();;

	return ret;
}

Ray Transform::operator()(const Ray& r) const
{
	Eigen::Vector3d o3 = (*this)(r.o);
	Eigen::Vector3d d3 = (*this)(r.d);

	return Ray(o3, d3, r.tMax, r.time, r.medium);
}

Bounds3f Transform::operator()(const Bounds3f& b) const
{

	const Transform& M = *this;
	Bounds3f ret(M(Eigen::Vector3d(b.pMin.x(), b.pMin.y(), b.pMin.z())));
	ret = Union(ret, M(Eigen::Vector3d(b.pMax.x(), b.pMin.y(), b.pMin.z())));
	ret = Union(ret, M(Eigen::Vector3d(b.pMin.x(), b.pMax.y(), b.pMin.z())));
	ret = Union(ret, M(Eigen::Vector3d(b.pMin.x(), b.pMin.y(), b.pMax.z())));
	ret = Union(ret, M(Eigen::Vector3d(b.pMin.x(), b.pMax.y(), b.pMax.z())));
	ret = Union(ret, M(Eigen::Vector3d(b.pMax.x(), b.pMax.y(), b.pMin.z())));
	ret = Union(ret, M(Eigen::Vector3d(b.pMax.x(), b.pMin.y(), b.pMax.z())));
	ret = Union(ret, M(Eigen::Vector3d(b.pMax.x(), b.pMax.y(), b.pMax.z())));
	return ret;


}

Transform Transform::operator*(const Transform& t2) const {
	return Transform(m * t2.m, mInv * t2.mInv);
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
	// TODO: p120
	return ret;
}
