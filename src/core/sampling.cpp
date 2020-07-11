#include "core/sampling.h"
#include "core/rng.h"
#include "utlis/utlis.h"


Distribution1D::Distribution1D(const Float* f, int n) : 
	func(f, f + n), cdf(n + 1)
{
	cdf[0] = 0;
	for (int i = 1; i < n + 1; ++i)
		cdf[i] = cdf[i - 1] + func[i - 1] / n;

	funcInt = cdf[n];
	if (funcInt == 0) {
		for (int i = 1; i < n + 1; ++i)
			cdf[i] = Float(i) / Float(n);
	}
	else {
		for (int i = 1; i < n + 1; ++i)
			cdf[i] /= funcInt;
	}
}

Float Distribution1D::SampleContinuous(Float u, Float* pdf, int* off) const
{
	int offset = FindInterval(cdf.size(),
		[&](int index) { return cdf[index] <= u; });

	if (off) *off 
		= offset;

	Float du = u - cdf[offset];
	if ((cdf[offset + 1] - cdf[offset]) > 0)
		du /= (cdf[offset + 1] - cdf[offset]);

	if (pdf) *pdf = func[offset] / funcInt;

	return (offset + du) / Count();
}

int Distribution1D::SampleDiscrete(Float u, Float* pdf = nullptr,
	Float* uRemapped = nullptr) const 
{
	int offset = FindInterval(cdf.size(),
		[&](int index) { return cdf[index] <= u; });

	if (pdf) 
		*pdf = func[offset] / (funcInt * Count());
	
	if (uRemapped)
		*uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
	return offset;
}

Float Distribution1D::DiscretePDF(int index) const {
	return func[index] / (funcInt * Count());
}

Distribution2D::Distribution2D(const Float* func, int nu, int nv)
{
	for (int v = 0; v < nv; ++v)
	{
		pConditionalV.emplace_back(new Distribution1D(&func[v * nu], nu));
	}
	std::vector<Float> marginalFunc;
	for (int v = 0; v < nv; ++v)
		marginalFunc.push_back(pConditionalV[v]->funcInt);
	pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));

}

Point2f Distribution2D::SampleContinuous(const Point2f& u, Float* pdf) const
{
	Float pdfs[2];
	int v;
	Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
	Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
	*pdf = pdfs[0] * pdfs[1];
	return Point2f(d0, d1);
}

Float Distribution2D::Pdf(const Point2f& p) const
{
	int iu = Clamp(int(p[0] * pConditionalV[0]->Count()),
		0, pConditionalV[0]->Count() - 1);
	int iv = Clamp(int(p[1] * pMarginal->Count()),
		0, pMarginal->Count() - 1);
	return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
}


Point2f RejectionSampleDisk(RNG& rng) 
{
	Point2f p;
	do {
		p.x() = 1 - 2 * rng.UniformFloat();
		p.y() = 1 - 2 * rng.UniformFloat();
	} while (p.x() * p.x() + p.y() * p.y() > 1);
	return p;
}

Vector3f UniformSampleHemisphere(const Point2f& u) 
{
	Float z = u[0];
	Float r = std::sqrt(std::max((Float)0, (Float)1. - z * z));
	Float phi = 2 * Pi * u[1];
	return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

Float UniformHemispherePdf() 
{
	return Inv2Pi;
}

Point2f UniformSampleDisk(const Point2f& u) 
{
	Float r = std::sqrt(u[0]);
	Float theta = 2 * Pi * u[1];
	return Point2f(r * std::cos(theta), r * std::sin(theta));
}

Point2f ConcentricSampleDisk(const Point2f& u) 
{
	Point2f uOffset = 2.0 * u - Vector2f(1, 1);

	if (uOffset.x() == 0 && uOffset.y() == 0)
		return Point2f(0, 0);

	Float theta, r;
	if (std::abs(uOffset.x()) > std::abs(uOffset.y())) 
	{
		r = uOffset.x();
		theta = PiOver4 * (uOffset.y() / uOffset.x());
	}
	else {
		r = uOffset.y();
		theta = PiOver2 - PiOver4 * (uOffset.x() / uOffset.y());
	}
	return r * Point2f(std::cos(theta), std::sin(theta));
}

Float UniformConePdf(Float cosThetaMax)
{
	return 1 / (2 * Pi * (1 - cosThetaMax));
}

Vector3f UniformSampleCone(const Point2f& u, Float cosThetaMax) 
{
	Float cosTheta = ((Float)1 - u[0]) + u[0] * cosThetaMax;
	Float sinTheta = std::sqrt((Float)1 - cosTheta * cosTheta);
	Float phi = u[1] * 2 * Pi;
	return Vector3f(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta,
		cosTheta);
}

Point2f UniformSampleTriangle(const Point2f& u)
{
	Float su0 = std::sqrt(u[0]);
	return Point2f(1 - su0, u[1] * su0);
}