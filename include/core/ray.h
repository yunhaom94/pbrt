#pragma once
#include <Eigen/Core>
#include <limits>

#include "core/medium.h"

#define inf std::numeric_limits<double>::infinity()

class Ray
{
public:
	// r(t) = o + td
	Eigen::Vector3d o;
	Eigen::Vector3d d;
	mutable double tMax;
	double time;
	const Medium *medium;

public:
	Ray() : tMax(inf), time(0), medium(nullptr) {}
	Ray(const Eigen::Vector3d& o, const Eigen::Vector3d& d, double tMax = inf,
		double time = 0, const Medium* medium = nullptr)
		: o(o), d(d), tMax(tMax), time(time), medium(medium) {}

	~Ray();

	// overload () operator
	Eigen::Vector3d operator()(double t) const;

private:

};

