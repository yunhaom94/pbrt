#include "core/ray.h"

Point3f Ray::operator()(Float t) const
{
	return o + d * t;
}

void RayDifferential::ScaleDifferentials(Float s)
{
	rxOrigin = o + (rxOrigin - o) * s;
	ryOrigin = o + (ryOrigin - o) * s;
	rxDirection = d + (rxDirection - d) * s;
	ryDirection = d + (ryDirection - d) * s;
}
