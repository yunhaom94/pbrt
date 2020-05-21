#include "cameras/camera.h"

Camera::Camera()
{
}

double Camera::GenerateRayDifferential(const CameraSample& sample, RayDifferential* rd) const
{
	
		double wt = GenerateRay(sample, rd);

		CameraSample sshift = sample;
		sshift.pFilm.x++;
		Ray rx;
		double wtx = GenerateRay(sshift, &rx);
		if (wtx == 0)
			return 0;
		rd->rxOrigin = rx.o;
		rd->rxDirection = rx.d;

		Ray ry;
		double wty = GenerateRay(sshift, &ry);
		if (wty == 0)
			return 0;
		rd->ryOrigin = ry.o;
		rd->ryDirection = ry.d;

		rd->hasDifferentials = true;
		return wt;
	
}

Camera::~Camera()
{
}
