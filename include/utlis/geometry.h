#pragma once
#include "core/pbrt.h"
#include "utlis/utlis.h"

// convert Spherical coordinate to Cartesian 
inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi)
{
	return Vector3f(sinTheta * std::cos(phi),
					sinTheta * std::sin(phi),
					cosTheta);
}

// // convert Spherical coordinate to Cartesian 
inline Eigen::Vector3d SphericalDirection(Float sinTheta, Float cosTheta, Float phi,
	const Vector3f& x, const Vector3f& y, const Vector3f& z)
{
	return sinTheta * std::cos(phi) * x + 
		   sinTheta * std::sin(phi) * y +
		   cosTheta * z;
}

inline Float SphericalTheta(const Vector3f& v)
{
	return std::acos(Clamp(v.z(), -1, 1));
}

inline Float SphericalPhi(const Vector3f& v) {
	Float p = std::atan2(v.y(), v.x());
	return (p < 0) ? (p + 2 * Pi) : p;
}


template <typename T>
int MaxDimension(const Vector3<T>& v) 
{
	return (v.x() > v.y()) ? ((v.x() > v.z()) ? 0 : 2) : ((v.y() > v.z()) ? 1 : 2);
}

template <typename T>
Vector3<T> Permute(const Vector3<T>& v, int x, int y, int z) {
	return Vector3<T>(v[x], v[y], v[z]);
}

template <typename T>
inline void CoordinateSystem(const Vector3<T>& v1, Vector3<T>* v2,
	Vector3<T>* v3) {
	if (std::abs(v1.x()) > std::abs(v1.y()))
		*v2 = Vector3<T>(-v1.z(), 0, v1.x()) / std::sqrt(v1.x() * v1.x() + v1.z() * v1.z());
	else
		*v2 = Vector3<T>(0, v1.z(), -v1.y()) / std::sqrt(v1.y() * v1.y() + v1.z() * v1.z());
	*v3 = v1.cross(*v2);
}