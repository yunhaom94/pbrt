#pragma once
#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/transport_mode.h"
#include "utlis/bsdf.h"

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

	Spectrum Evaluate(Float cosThetaI) const;

};

class FresnelDielectric : public Fresnel 
{
private:
	Float etaI, etaT;

public:
	FresnelDielectric(Float etaI, Float etaT) : etaI(etaI), etaT(etaT) { }

	Spectrum Evaluate(Float cosThetaI) const;

};

class FresnelNoOp : public Fresnel 
{
public:
	Spectrum Evaluate(Float) const { return Spectrum(1.0); }
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
	// (defused)
	virtual Spectrum f(const Vector3f& wo, const Vector3f& wi) const = 0;

	// with only one direction (specular)
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

class SpecularReflection : public BxDF 
{
private:
	const Spectrum R;
	const Fresnel* fresnel;

public:
	SpecularReflection(const Spectrum& R, Fresnel* fresnel)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)), R(R),
		fresnel(fresnel) { }

	// for an arbitrary pair of directions the delta function
	// returns no scattering (specular only).
	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;

	// specular
	Spectrum Sample_f(const Vector3f& wo,
		Vector3f* wi, const Point2f& sample, Float* pdf,
		BxDFType* sampledType) const;
};


class SpecularTransmission : public BxDF
{
private:
	const Spectrum T;
	const Float etaA; // index of refraction above the surface 
	const Float etaB; // index of refraction below the surface 
	const FresnelDielectric fresnel;
	const TransportMode mode;

public:
	SpecularTransmission(const Spectrum& T, Float etaA, Float etaB,
		TransportMode mode)
		: BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)), T(T), etaA(etaA),
		etaB(etaB), fresnel(etaA, etaB), mode(mode) { }

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;

	Spectrum Sample_f(const Vector3f& wo,
		Vector3f* wi, const Point2f& sample, Float* pdf,
		BxDFType* sampledType) const;
};

class FresnelSpecular : public BxDF 
{
public:
	const Spectrum R, T;
	const Float etaA, etaB;
	const FresnelDielectric fresnel;
	const TransportMode mode;

public:
	FresnelSpecular(const Spectrum& R, const Spectrum& T, Float etaA,
		Float etaB, TransportMode mode)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)),
		R(R), T(T), etaA(etaA), etaB(etaB), fresnel(etaA, etaB),
		mode(mode) { }

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;
};

// non-accurate reflection model
// scatters incident illumination equally in all directions.
// what we used in 418
class LambertianReflection : public BxDF 
{
private:
	const Spectrum R;

public:
	LambertianReflection(const Spectrum& R)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) { }

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;

	Spectrum rho(const Vector3f&, int, const Point2f*) const { return R; }
	Spectrum rho(int, const Point2f*, const Point2f*) const { return R; }

};


class OrenNayar : public BxDF
{
private:
	const Spectrum R;
	Float A, B;

public:
	OrenNayar(const Spectrum& R, Float sigma)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R)
	{
		sigma = Radians(sigma);
		Float sigma2 = sigma * sigma;
		A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
		B = 0.45f * sigma2 / (sigma2 + 0.09f);
	}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;


};


class MicrofacetReflection : public BxDF
{
private:
	const Spectrum R;
	const MicrofacetDistribution* distribution;
	const Fresnel* fresnel;

public:
	MicrofacetReflection(const Spectrum & R,
		MicrofacetDistribution * distribution, Fresnel * fresnel)
	: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), R(R),
	distribution(distribution), fresnel(fresnel) { }


	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;

};

class MicrofacetTransmission : public BxDF 
{
private:
	const Spectrum T;
	const MicrofacetDistribution* distribution;
	const Float etaA, etaB;
	const FresnelDielectric fresnel;
	const TransportMode mode;

public:
	MicrofacetTransmission(const Spectrum& T,
		MicrofacetDistribution* distribution, Float etaA, Float etaB,
		TransportMode mode)
		: BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_GLOSSY)),
		T(T), distribution(distribution), etaA(etaA), etaB(etaB),
		fresnel(etaA, etaB), mode(mode) { }

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;
};


// surface with glossy on top of diffuse
class FresnelBlend : public BxDF
{
private:
	const Spectrum Rd, Rs;
	MicrofacetDistribution* distribution;

public:
	FresnelBlend::FresnelBlend(const Spectrum& Rd, const Spectrum& Rs,
		MicrofacetDistribution* distribution)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
		Rd(Rd), Rs(Rs), distribution(distribution) { }

	Spectrum SchlickFresnel(Float cosTheta) const {
		auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
		return Rs + pow5(1 - cosTheta) * (Spectrum(1.) - Rs);
	}
	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;
};


struct FourierBSDFTable 
{
	Float eta;
	int mMax;
	int nChannels;
	int nMu;
	Float* mu;
	int* m;
	int* aOffset;
	Float* a;

	static bool Read(const std::string& filename, FourierBSDFTable* table);
	const Float* GetAk(int offsetI, int offsetO, int* mptr) const;

	bool GetWeightsAndOffset(Float cosTheta, int* offset, Float weights[4]) const;

};

class FourierBSDF : public BxDF
{

private:
	const FourierBSDFTable& bsdfTable;
	const TransportMode mode;

public:
	FourierBSDF(const FourierBSDFTable& bsdfTable, TransportMode mode)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)),
		bsdfTable(bsdfTable), mode(mode) { }

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const;
};