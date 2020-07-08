#pragma once
#pragma once
#include "core/light.h"
#include "core/texture.h"

// TODO: p726
class ProjectionLight : public Light 
{
	/*
private:
	Point2i resolution;
	std::unique_ptr<RGBSpectrum[]> texels = ReadImage(texname, &resolution);
	if (texels)
		projectionMap.reset(new MIPMap<RGBSpectrum>(resolution,
			texels.get()));

public:

	ProjectionLight(const Transform& LightToWorld, const MediumInterface& mediumInterface, const Spectrum& I, const std::string& texname, Float fov);
	*/
};