
#include "core/transform.h"



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

inline Ray Transform::operator()(const Ray& r) const
{

	Eigen::Vector4d o = this->m * Eigen::Vector4d(r.o.x(), r.o.y(), r.o.z(), 0);
	Eigen::Vector4d d = this->m * Eigen::Vector4d(r.d.x(), r.d.y(), r.d.z(), 0);

	Eigen::Vector3d o3 = o.head<3>();
	Eigen::Vector3d d3 = d.head<3>();

	return Ray(o3, d3, r.tMax, r.time, r.medium);
}

Transform Transform::operator*(const Transform& t2) const {
	return Transform(m * t2.m, mInv * t2.mInv);
}