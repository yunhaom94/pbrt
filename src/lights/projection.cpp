#include "lights/projection.h"

/*
ProjectionLight::ProjectionLight(const Transform& LightToWorld,
	const MediumInterface& mediumInterface, const Spectrum& I,
	const std::string& texname, Float fov)
	: Light((int)LightFlags::DeltaPosition, LightToWorld,
		mediumInterface),
	pLight(LightToWorld(Point3f(0, 0, 0))), I(I)
{
	Create ProjectionLight MIP map 726
		Initialize ProjectionLight projection matrix 726
		Compute cosine of cone surrounding projection directions 727
}
*/