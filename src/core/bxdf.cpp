#include "core/bxdf.h"
#include "core/microfacet.h"

// Fresnel reflection formula for dielectric
// materialsand unpolarized light.
// Used to compute amount of light reflected
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

Spectrum SpecularReflection::f(const Vector3f& wo, const Vector3f& wi) const
{
	return Spectrum(0);
}

Spectrum SpecularReflection::Sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, Float* pdf, BxDFType* sampledType) const
{

	*wi = Vector3f(-wo.x(), -wo.y(), wo.z());

	*pdf = 1;
	return fresnel->Evaluate(CosTheta(*wi)) * R / AbsCosTheta(*wi);
}

Spectrum SpecularTransmission::f(const Vector3f& wo, const Vector3f& wi) const
{
	return Spectrum(0);
}

Spectrum SpecularTransmission::Sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, Float* pdf, BxDFType* sampledType) const
{
	bool entering = CosTheta(wo) > 0;
	Float etaI = entering ? etaA : etaB;
	Float etaT = entering ? etaB : etaA;

	if (!Refract(wo, Faceforward(Normal3f(0, 0, 1), wo), etaI / etaT, wi))
		return 0;

	* pdf = 1;
	Spectrum ft = T * (Spectrum(1.) - fresnel.Evaluate(CosTheta(*wi)));
	// TODO: p961 Account for non - symmetry with transmission to different medium 961
	return ft / AbsCosTheta(*wi);
}

Spectrum FresnelSpecular::f(const Vector3f& wo, const Vector3f& wi) const
{
	return Spectrum(0.0);
}

Spectrum FresnelConductor::Evaluate(Float cosThetaI) const
{
	return FrConductor(std::abs(cosThetaI), etaI, etaT, k);
}

Spectrum FresnelDielectric::Evaluate(Float cosThetaI) const
{
	return FrDielectric(cosThetaI, etaI, etaT);
}

Spectrum LambertianReflection::f(const Vector3f& wo, const Vector3f& wi) const
{
	return R * InvPi;
}

Spectrum OrenNayar::f(const Vector3f& wo, const Vector3f& wi) const
{
	Float sinThetaI = SinTheta(wi);
	Float sinThetaO = SinTheta(wo);

	Float maxCos = 0;
	if (sinThetaI > 1e-4 && sinThetaO > 1e-4) {
		Float sinPhiI = SinPhi(wi), cosPhiI = CosPhi(wi);
		Float sinPhiO = SinPhi(wo), cosPhiO = CosPhi(wo);
		Float dCos = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
		maxCos = std::max((Float)0, dCos);
	}

	Float sinAlpha, tanBeta;
	if (AbsCosTheta(wi) > AbsCosTheta(wo)) {
		sinAlpha = sinThetaO;
		tanBeta = sinThetaI / AbsCosTheta(wi);
	}
	else {
		sinAlpha = sinThetaI;
		tanBeta = sinThetaO / AbsCosTheta(wo);
	}

	return R * InvPi * (A + B * maxCos * sinAlpha * tanBeta);
}

Spectrum MicrofacetReflection::f(const Vector3f& wo,
	const Vector3f& wi) const {
	Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
	Vector3f wh = wi + wo;
	
	if (cosThetaI == 0 || cosThetaO == 0) 
		return Spectrum(0);

	if (wh.x() == 0 && wh.y() == 0 && wh.z() == 0)
		return Spectrum(0);

	wh = wh.normalized();
	Spectrum F = fresnel->Evaluate(wi.dot(wh));
	return R * distribution->D(wh) * distribution->G(wo, wi) * F /
		(4 * cosThetaI * cosThetaO);
}

Spectrum MicrofacetTransmission::f(const Vector3f& wo, const Vector3f& wi) const
{
	if (SameHemisphere(wo, wi)) return 0;  // transmission only

	Float cosThetaO = CosTheta(wo);
	Float cosThetaI = CosTheta(wi);
	if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0);

	// Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
	Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
	Vector3f wh = (wo + wi * eta).normalized();
	if (wh.z() < 0) wh = -wh;

	Spectrum F = fresnel.Evaluate(wo.dot(wh));

	Float sqrtDenom = wo.dot(wh) + eta * wi.dot(wh);
	Float factor = (mode == TransportMode::Radiance) ? (1 / eta) : 1;

	return (Spectrum(1.f) - F) * T *
		std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta *
			std::abs(wi.dot(wh)) * std::abs(wo.dot(wh)) * factor * factor /
			(cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
}

Spectrum FresnelBlend::f(const Vector3f& wo, const Vector3f& wi) const {
	auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
	Spectrum diffuse = (28.f / (23.f * Pi)) * Rd *
		(Spectrum(1.f) - Rs) *
		(1 - pow5(1 - .5f * AbsCosTheta(wi))) *
		(1 - pow5(1 - .5f * AbsCosTheta(wo)));
	Vector3f wh = wi + wo;
	
	if (wh.x() == 0 && wh.y() == 0 && wh.z() == 0) 
		return Spectrum(0);
	
	wh = wh.normalized();
	Spectrum specular = distribution->D(wh) /
		(4 * std::abs(wi.dot(wh)) *
			std::max(AbsCosTheta(wi), AbsCosTheta(wo))) *
		SchlickFresnel(wi.dot(wh));
	return diffuse + specular;
}

bool FourierBSDFTable::Read(const std::string& filename, FourierBSDFTable* table)
{
	// TODO:
	return false;
}

Float Fourier(const Float* a, int m, double cosPhi) {
	double value = 0.0;
	double cosKMinusOnePhi = cosPhi;
	double cosKPhi = 1;
	for (int k = 0; k < m; ++k) {
		value += a[k] * cosKPhi;
		double cosKPlusOnePhi = 2 * cosPhi * cosKPhi - cosKMinusOnePhi;
		cosKMinusOnePhi = cosKPhi;
		cosKPhi = cosKPlusOnePhi;
	}
	return value;
}
const Float* FourierBSDFTable::GetAk(int offsetI, int offsetO, int* mptr) const
{
	*mptr = m[offsetO * nMu + offsetI];
	return a + aOffset[offsetO * nMu + offsetI];
}



bool FourierBSDFTable::GetWeightsAndOffset(Float cosTheta, int* offset,
	Float weights[4]) const {
	return CatmullRomWeights(nMu, mu, cosTheta, offset, weights);
}

Spectrum FourierBSDF::f(const Vector3f& wo, const Vector3f& wi) const
{
	Float muI = CosTheta(-wi), muO = CosTheta(wo);
	Float cosPhi = CosDPhi(-wi, wo);

	int offsetI, offsetO;
	Float weightsI[4], weightsO[4];
	if (!bsdfTable.GetWeightsAndOffset(muI, &offsetI, weightsI) ||
		!bsdfTable.GetWeightsAndOffset(muO, &offsetO, weightsO))
		return Spectrum(0);

	Float* ak = ALLOCA(Float, bsdfTable.mMax * bsdfTable.nChannels);
	memset(ak, 0, bsdfTable.mMax * bsdfTable.nChannels * sizeof(Float));

	int mMax = 0;
	for (int b = 0; b < 4; ++b) {
		for (int a = 0; a < 4; ++a) {
			Float weight = weightsI[a] * weightsO[b];
			if (weight != 0) {
				int m;
				const Float* ap = bsdfTable.GetAk(offsetI + a, offsetO + b, &m);
				mMax = std::max(mMax, m);
				for (int c = 0; c < bsdfTable.nChannels; ++c)
					for (int k = 0; k < m; ++k)
						ak[c * bsdfTable.mMax + k] += weight * ap[c * m + k];
			}
		}
	}

	Float Y = std::max((Float)0, Fourier(ak, mMax, cosPhi));
	Float scale = muI != 0 ? (1 / std::abs(muI)) : (Float)0;
	// TODO: p961 Update scale to account for adjoint light transport 961
		if (bsdfTable.nChannels == 1)
			return Spectrum(Y * scale);
		else {

			Float R = Fourier(ak + 1 * bsdfTable.mMax, mMax, cosPhi);
			Float B = Fourier(ak + 2 * bsdfTable.mMax, mMax, cosPhi);
			Float G = 1.39829f * Y - 0.100913f * B - 0.297375f * R;
			Float rgb[3] = { R * scale, G * scale, B * scale };
			return Spectrum::FromRGB(rgb).Clamp();
		}
}

Spectrum BSDF::f(const Vector3f& woW, const Vector3f& wiW,
	BxDFType flags) const {
	Vector3f wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
	bool reflect = wiW.dot(ng) * woW.dot(ng) > 0;
	Spectrum f(0.f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags) &&
			((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
				(!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
			f += bxdfs[i]->f(wo, wi);
	return f;
}

Spectrum BSDF::Sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& u, Float* pdf, BxDFType* sampledType) const
{
	// TODO:
	return Spectrum();
}

Spectrum BSDF::rho(int nSamples, const Point2f* samples1, const Point2f* samples2, BxDFType flags) const
{
	Spectrum ret(0.f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags))
			ret += bxdfs[i]->rho(nSamples, samples1, samples2);
	return ret;
}

Spectrum BSDF::rho(const Vector3f& wo, int nSamples, const Point2f* samples, BxDFType flags) const
{
	Vector3f wo = WorldToLocal(wo);
	Spectrum ret(0.f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags))
			ret += bxdfs[i]->rho(wo, nSamples, samples);
	return ret;
}
