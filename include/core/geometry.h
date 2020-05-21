#pragma once
#include "core/pbrt.h"

// using inline so that I can define it in .h
inline Eigen::Vector3d SphericalDirection(double sinTheta,
	double cosTheta, double phi) {
	return Eigen::Vector3d(sinTheta * std::cos(phi),
		sinTheta * std::sin(phi),
		cosTheta);
}

inline Eigen::Vector3d SphericalDirection(double sinTheta, double cosTheta, double phi,
	const Eigen::Vector3d& x, const Eigen::Vector3d& y,
	const Eigen::Vector3d& z) {
	return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y +
		cosTheta * z;
}

inline double SphericalTheta(const Eigen::Vector3d& v) {
	return std::acos(Clamp<double, double, double>(v.z(), -1, 1));
}

inline double SphericalPhi(const Eigen::Vector3d& v) {
	double p = std::atan2(v.y(), v.x());
	return (p < 0) ? (p + 2 * Pi) : p;
}