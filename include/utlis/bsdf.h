#pragma once
#include "core/pbrt.h"
#include "utlis/utlis.h"

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

inline Vector3f Reflect(const Vector3f& wo, const Vector3f& n) 
{
	return -wo + 2 * wo.dot(n) * n;
}

inline Float CosPhi(const Vector3f& w) {
	Float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 1 : Clamp(w.x() / sinTheta, -1, 1);
}
inline Float SinPhi(const Vector3f& w) {
	Float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 0 : Clamp(w.y() / sinTheta, -1, 1);
}

inline Float CosDPhi(const Vector3f& wa, const Vector3f& wb) {
	return Clamp((wa.x() * wb.x() + wa.y() * wb.y()) /
		std::sqrt((wa.x() * wa.x() + wa.y() * wa.y()) *
			(wb.x() * wb.x() + wb.y() * wb.y())), -1, 1);
}

inline Float Cos2Phi(const Vector3f& w) {
	return CosPhi(w) * CosPhi(w);
}
inline Float Sin2Phi(const Vector3f& w) {
	return SinPhi(w) * SinPhi(w);
}

inline bool Refract(const Vector3f& wi, const Normal3f& n, Float eta,
	Vector3f* wt)
{
	// Snell¡¯s law
	Float cosThetaI = n.dot(wi);
	Float sin2ThetaI = std::max(0.0, 1 - cosThetaI * cosThetaI);
	Float sin2ThetaT = eta * eta * sin2ThetaI;
	
	// total internal reflection
	if (sin2ThetaT >= 1) return
		false;
		
	Float cosThetaT = std::sqrt(1 - sin2ThetaT);

	*wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
	return true;
}

inline bool SameHemisphere(const Vector3f& w, const Vector3f& wp) {
	return w.z() * wp.z() > 0;
}

inline bool SameHemisphere(const Vector3f& w, const Normal3f& wp) {
	return w.z() * wp.z() > 0;
}