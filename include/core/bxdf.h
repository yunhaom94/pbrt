#pragma once
#include "core/pbrt.h"
#include "core/spectrum.h"

enum BxDFType {
	BSDF_REFLECTION = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
	BSDF_DIFFUSE = 1 << 2,
	BSDF_GLOSSY = 1 << 3,
	BSDF_SPECULAR = 1 << 4,
	BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR |
	BSDF_REFLECTION | BSDF_TRANSMISSION,
};

Float FrDielectric(Float cosThetaI, Float etaI, Float etaT);

Spectrum FrConductor(Float cosThetaI, const Spectrum& etaI,
	const Spectrum& etaT, const Spectrum& k);

class Fresnel 
{
public:
	virtual Spectrum Evaluate(Float cosI) const = 0;
};

class FresnelConductor : public Fresnel 
{
private:
	Spectrum etaI, etaT, k;

public:
	FresnelConductor(const Spectrum& etaI, const Spectrum& etaT,
		const Spectrum& k) : etaI(etaI), etaT(etaT), k(k) { }

};

class FresnelDielectric : public Fresnel 
{
private:
	Float etaI, etaT;

public:
	FresnelDielectric(Float etaI, Float etaT) : etaI(etaI), etaT(etaT) { }

};

class FresnelNoOp : public Fresnel 
{
public:
	Spectrum Evaluate(Float) const { return Spectrum(1); }
};

// describe a type of light behavior on a point on a surface
class BxDF
{
public:
	const BxDFType type;


public:
	BxDF(BxDFType type) : type(type) { }
	~BxDF() {};

	bool MatchesFlags(BxDFType t) const 
	{
		return (type & t) == type;
	}

	// returns the value of the distribution function for the given pair of directions.
	virtual Spectrum f(const Vector3f& wo, const Vector3f& wi) const = 0;

	// with only one direction
	virtual Spectrum Sample_f(const Vector3f& wo, Vector3f* wi,
		const Point2f& sample, Float* pdf,
		BxDFType* sampledType = nullptr) const;

	// hemispherical-directional reflectance function
	virtual Spectrum rho(const Vector3f& wo, int nSamples,
		const Point2f* samples) const;

	// compute hemispherical-hemispherical reflectance of a surface when 
	// no w0 given
	virtual Spectrum rho(int nSamples, const Point2f* samples1,
		const Point2f* samples2) const;

private:

};

// take a given BxDF and scale its contribution with a Spectrum value.
class ScaledBxDF : public BxDF 
{

private:
	BxDF* bxdf;
	Spectrum scale;

public:
	ScaledBxDF(BxDF* bxdf, const Spectrum& scale)
		: BxDF(BxDFType(bxdf->type)), bxdf(bxdf), scale(scale) {}


	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;

	// with only one direction
	virtual Spectrum Sample_f(const Vector3f& wo, Vector3f* wi,
		const Point2f& sample, Float* pdf,
		BxDFType* sampledType = nullptr) const;

	virtual Spectrum rho(const Vector3f& wo, int nSamples,
		const Point2f* samples) const;

	virtual Spectrum rho(int nSamples, const Point2f* samples1,
		const Point2f* samples2) const;

};