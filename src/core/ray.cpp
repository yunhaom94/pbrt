#include "core/ray.h"

Eigen::Vector3d Ray::operator()(double t) const
{
	return o + d * t;
}