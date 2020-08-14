#pragma once

#include "core/pbrt.h"
#include "utlis/bsdf.h"

// abstraction for more advance
// Microfacet distributions functions
class MicrofacetDistribution
{
protected:
	const bool sampleVisibleArea;

public:

	~MicrofacetDistribution() {}

	// return differential area of microfacets 
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

	virtual Vector3f Sample_wh(const Vector3f& wo,
		const Point2f& u) const = 0;

	Float Pdf(const Vector3f& wo, const Vector3f& wh) const;

protected:
	MicrofacetDistribution(bool sampleVisibleArea);

};

class BeckmannDistribution : public MicrofacetDistribution 
{
private:
	const Float alphax, alphay;

public:
	BeckmannDistribution(Float alphax, Float alphay, bool samplevis = true)
		: MicrofacetDistribution(samplevis), alphax(alphax), alphay(alphay) {}

	Float D(const Vector3f& wh) const;

	Float Lambda(const Vector3f& w) const;

	Vector3f Sample_wh(const Vector3f& wo,
		const Point2f& u) const;
	
};

class TrowbridgeReitzDistribution : public MicrofacetDistribution 
{

private:
	const Float alphax, alphay;

public:
	TrowbridgeReitzDistribution(Float alphax, Float alphay, bool samplevis = true)
		: MicrofacetDistribution(samplevis), alphax(alphax), alphay(alphay) {}

	static inline Float RoughnessToAlpha(Float roughness) { return 0; }// TODO:

	Float D(const Vector3f& wh) const;
	Float Lambda(const Vector3f& w) const;

	virtual Vector3f Sample_wh(const Vector3f& wo,
		const Point2f& u) const;
};