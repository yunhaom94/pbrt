#pragma once
#include "core/pbrt.h"
#include "utlis/bsdf.h"

//Given a direction vector ¦Ø in this coordinate system, it is easy to compute quantities like
// the cosine of the angle that it forms with the normal direction
inline double CosTheta(const Eigen::Vector3d& w) { return w.z(); }
// cos^2theta
inline double Cos2Theta(const Eigen::Vector3d& w) { return w.z() * w.z(); }
inline double AbsCosTheta(const Eigen::Vector3d& w) { return std::abs(w.z()); }

// sin^2theta
inline double Sin2Theta(const Eigen::Vector3d& w)
{
	return std::max(0.0, 1.0 - Cos2Theta(w));
}
inline double SinTheta(const Eigen::Vector3d& w)
{
	return std::sqrt(Sin2Theta(w));
}

inline double TanTheta(const Eigen::Vector3d& w)
{
	return SinTheta(w) / CosTheta(w);
}
inline double Tan2Theta(const Eigen::Vector3d& w) 
{
	return Sin2Theta(w) / Cos2Theta(w);
}