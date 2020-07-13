#pragma once

#include "core/pbrt.h"

struct Distribution1D 
{
	std::vector<Float> func, cdf;
	Float funcInt;

	Distribution1D(const Float* f, int n);
	~Distribution1D() {}

	int Count() const { return func.size(); }

	Float SampleContinuous(Float u, Float* pdf, int* off = nullptr) const;

	int SampleDiscrete(Float u, Float* pdf, Float* uRemapped) const;

	Float DiscretePDF(int index) const;

};

class Distribution2D
{
private:
	std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
	std::unique_ptr<Distribution1D> pMarginal;

public:
	Distribution2D(const Float* func, int nu, int nv);

	~Distribution2D() {}

	Point2f SampleContinuous(const Point2f& u, Float* pdf) const;

	Float Pdf(const Point2f& p) const;


};


Point2f RejectionSampleDisk(RNG& rng);

Vector3f UniformSampleHemisphere(const Point2f& u);

Float UniformHemispherePdf();

Point2f UniformSampleDisk(const Point2f& u);

Point2f ConcentricSampleDisk(const Point2f& u);

inline Vector3f CosineSampleHemisphere(const Point2f& u)
{
	Point2f d = ConcentricSampleDisk(u);
	Float z = std::sqrt(std::max((Float)0, 1 - d.x() * d.x() - d.y() * d.y()));
	return Vector3f(d.x(), d.y(), z);
}

inline Float CosineHemispherePdf(Float cosTheta) 
{
	return cosTheta * InvPi;
}

Float UniformConePdf(Float cosThetaMax);

Vector3f UniformSampleCone(const Point2f& u, Float cosThetaMax);

Point2f UniformSampleTriangle(const Point2f& u);

inline Float BalanceHeuristic(int nf, Float fPdf, int ng, Float gPdf) 
{
	return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}

inline Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf)
{
	Float f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}


Vector3f UniformSampleSphere(const Point2f& u);

Float UniformSpherePdf();