#pragma once

#include "core/pbrt.h"
#include "utlis/bsdf.h"

// faces distributions
class MicrofacetDistribution
{
public:
	MicrofacetDistribution() {}
	~MicrofacetDistribution() {}


	// the differential area of microfacets 
	// oriented with the given normal vector
	virtual Float D(const Vector3f& wh) const = 0;

	// shadow masking of microfacets
	virtual Float Lambda(const Vector3f& w) const = 0;

	Float G1(const Vector3f& w) const
	{
		return 1 / (1 + Lambda(w));
	}

	Float G(const Vector3f& wo, const Vector3f& wi) const
	{
		return 1 / (1 + Lambda(wo) + Lambda(wi));
	}

protected:

};

class BeckmannDistribution : public MicrofacetDistribution 
{
private:
	const Float alphax, alphay;

public:
	BeckmannDistribution(Float alphax, Float alphay)
		: alphax(alphax), alphay(alphay) {}

	Float D(const Vector3f& wh) const;

	Float Lambda(const Vector3f& w) const;

	
};

class TrowbridgeReitzDistribution : public MicrofacetDistribution 
{

private:
	const Float alphax, alphay;

public:
	TrowbridgeReitzDistribution(Float alphax, Float alphay)
		: alphax(alphax), alphay(alphay) {}

	static inline Float RoughnessToAlpha(Float roughness) { return 0; }// TODO:

	Float D(const Vector3f& wh) const;
	Float Lambda(const Vector3f& w) const;
};