#include "core/bxdf.h"

// Fresnel reflection formula for dielectric
// materialsand unpolarized light.
Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) 
{
	cosThetaI = Clamp(cosThetaI, -1, 1);

	bool entering = cosThetaI > 0.f;
	if (!entering) {
		std::swap(etaI, etaT);
		cosThetaI = std::abs(cosThetaI);
	}

	Float sinThetaI = std::sqrt(std::max((Float)0,
		1 - cosThetaI * cosThetaI));
	Float sinThetaT = etaI / etaT * sinThetaI;
	
	if (sinThetaT >= 1)
		return 1;
	
	Float cosThetaT = std::sqrt(std::max((Float)0,
			1 - sinThetaT * sinThetaT));

	Float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
		((etaT * cosThetaI) + (etaI * cosThetaT));
	Float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
		((etaI * cosThetaI) + (etaT * cosThetaT));
	return (Rparl * Rparl + Rperp * Rperp) / 2;
}

Spectrum FrConductor(Float cosThetaI, const Spectrum& etaI, const Spectrum& etaT, const Spectrum& k)
{
	// TODO:
	return Spectrum();
}

Spectrum ScaledBxDF::f(const Vector3f& wo, const Vector3f& wi) const 
{
	return scale * bxdf->f(wo, wi);
}

Spectrum ScaledBxDF::Sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, Float* pdf, BxDFType* sampledType) const
{
	return scale * bxdf->Sample_f(wo, wi, sample, pdf, sampledType);
}

Spectrum ScaledBxDF::rho(const Vector3f& wo, int nSamples, const Point2f* samples) const
{
	return scale * bxdf->rho(wo, nSamples, samples);
}

Spectrum ScaledBxDF::rho(int nSamples, const Point2f* samples1, const Point2f* samples2) const
{
	return scale * bxdf->rho(nSamples, samples1, samples2);
}

