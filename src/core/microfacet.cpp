#include "core/microfacet.h"

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