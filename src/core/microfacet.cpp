#include "core/microfacet.h"
#include "utlis/geometry.h"

MicrofacetDistribution::MicrofacetDistribution(bool sampleVisibleArea)
	: sampleVisibleArea(sampleVisibleArea) { }

Float MicrofacetDistribution::Pdf(const Vector3f& wo,
	const Vector3f& wh) const 
{
	if (sampleVisibleArea)
		return D(wh) * G1(wo) * std::abs(wo.dot(wh)) / AbsCosTheta(wo);
	else
		return D(wh) * AbsCosTheta(wh);
}

Float BeckmannDistribution::D(const Vector3f& wh) const
{
	Float tan2Theta = Tan2Theta(wh);

	if (std::isinf(tan2Theta)) 
		return 0.;

	Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
	return std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) +
		Sin2Phi(wh) / (alphay * alphay))) /
		(Pi * alphax * alphay * cos4Theta);
}

Float BeckmannDistribution::Lambda(const Vector3f& w) const {
	Float absTanTheta = std::abs(TanTheta(w));

	if (std::isinf(absTanTheta)) 
		return 0.;

	Float alpha = std::sqrt(Cos2Phi(w) * alphax * alphax +
		Sin2Phi(w) * alphay * alphay);

	Float a = 1 / (alpha * absTanTheta);
	if (a >= 1.6f)
		return 0;
	return (1 - 1.259f * a + 0.396f * a * a) /
		(3.535f * a + 2.181f * a * a);
}

Vector3f BeckmannDistribution::Sample_wh(const Vector3f& wo, const Point2f& u) const
{
	if (!sampleVisibleArea)
	{
		Float tan2Theta, phi;
		if (alphax == alphay)
		{
			Float logSample = std::log(u[0]);
			if (std::isinf(logSample)) 
				logSample = 0;

			tan2Theta = -alphax * alphax * logSample;
			phi = u[1] * 2 * Pi;
		}
		else
		{
			Float logSample = std::log(1 - u[0]);

			phi = std::atan(alphay / alphax *
				std::tan(2 * Pi * u[1] + 0.5f * Pi));
			if (u[1] > 0.5f)
				phi += Pi;

			Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
			Float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
			tan2Theta = -logSample /
				(cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
		}
		Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
		Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
		Vector3f wh = SphericalDirection(sinTheta, cosTheta, phi);
		if (!SameHemisphere(wo, wh)) wh = -wh;
			return wh;
	}
	else 
	{
		// Sample visible area of normals for Beckmann distribution
		Vector3f wh;
		bool flip = wo.z() < 0;
		// TODO:
		//wh = BeckmannSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
		if (flip) 
			wh = -wh;
		return wh;
	}
}


Float TrowbridgeReitzDistribution::D(const Vector3f& wh) const
{
	Float tan2Theta = Tan2Theta(wh);
	if (std::isinf(tan2Theta)) 
		return 0.;
	const Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
	Float e = (Cos2Phi(wh) / (alphax * alphax) +
		Sin2Phi(wh) / (alphay * alphay)) * tan2Theta;
	return 1 / (Pi * alphax * alphay * cos4Theta * (1 + e) * (1 + e));
}


Float TrowbridgeReitzDistribution::Lambda(const Vector3f & w) const 
{
	Float absTanTheta = std::abs(TanTheta(w));
	if (std::isinf(absTanTheta)) 
		return 0.;
	
	Float alpha = std::sqrt(Cos2Phi(w) * alphax * alphax +
		Sin2Phi(w) * alphay * alphay);
		
	Float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
	return (-1 + std::sqrt(1.f + alpha2Tan2Theta)) / 2;
}

