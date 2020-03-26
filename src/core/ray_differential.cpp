#include "core\ray_differential.h"

void RayDifferential::ScaleDifferentials(double s)
{
	rxOrigin = o + (rxOrigin - o) * s;
	ryOrigin = o + (ryOrigin - o) * s;
	rxDirection = d + (rxDirection - d) * s;
	ryDirection = d + (ryDirection - d) * s;
}
