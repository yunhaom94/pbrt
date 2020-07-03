#include "core/interaction.h"
#include "core/medium_interface.h"
#include "core/ray.h"
#include "core/spectrum.h"
#include "core/shape.h"
#include "core/memory.h"
#include "core/primitive.h"
#include "utlis/utlis.h"

Interaction::Interaction(const Point3f& p, 
	const Normal3f& n, 
	const Vector3f& pError,
	const Vector3f& wo, 
	Float time,
	const MediumInterface mediumInterface) :
	p(p), 
	time(time), 
	pError(pError), 
	wo(wo),
	n(n),
	mediumInterface(mediumInterface) { }


bool Interaction::IsSurfaceInteraction() const
{
	return surfaceInteraction;
}

SurfaceInteraction::SurfaceInteraction(const Point3f& p,
	const Vector3f& pError, 
	const Point2f& uv, 
	const Vector3f& wo,
	const Vector3f& dpdu, 
	const Vector3f& dpdv,
	const Normal3f& dndu, 
	const Normal3f& dndv,
	Float time, 
	const Shape* shape)
	: Interaction(p, 
		dpdu.cross(dpdv).normalized(),
		pError,
		wo,
		time, 
		MediumInterface()), // TODO: this, figure out why book put nullptr
	uv(uv), 
	dpdu(dpdu), 
	dpdv(dpdv), 
	dndu(dndu),
	dndv(dndv),
	shape(shape) 
{
	surfaceInteraction = true;
	shading.n = n;
	shading.dpdu = dpdu;
	shading.dpdv = dpdv;
	shading.dndu = dndu;
	shading.dndv = dndv;

	// sometime normals needs to be reversed (specified in the input)
	if (shape && (shape->reverseOrientation ^ shape->transformSwapsHandedness)) 
	{
		n *= -1;
		shading.n *= -1;
	}
}



void SurfaceInteraction::SetShadingGeometry(const Vector3f& dpdus,
	const Vector3f& dpdvs,
	const Normal3f& dndus, 
	const Normal3f& dndvs, 
	bool orientationIsAuthoritative)
{
	shading.n = dpdus.cross(dpdvs).normalized();
	if (shape && (shape->reverseOrientation ^
		shape->transformSwapsHandedness))
		shading.n = -shading.n;

	// reverse the direction of normals if needed
	// so that they face the same direction 
	if (orientationIsAuthoritative)
		n = Faceforward(n, shading.n);
	else
		shading.n = Faceforward(shading.n, n);

	shading.dpdu = dpdus;
	shading.dpdv = dpdvs;
	shading.dndu = dndus;
	shading.dndv = dndvs;

}

void SurfaceInteraction::ComputeDifferentials(Ray r)
{
	// TODO: p601
}

void SurfaceInteraction::ComputeScatteringFunctions(const RayDifferential& ray, MemoryArena& arena, bool allowMultipleLobes, TransportMode mode)
{
	ComputeDifferentials(ray);
	primitive->ComputeScatteringFunctions(this, arena, mode,
		allowMultipleLobes);
}



Spectrum SurfaceInteraction::Le(Vector3f wo)
{
	// TODO:
	return Spectrum();
}
