#pragma once

#include "core/pbrt.h"
#include "core/medium_interface.h"

struct Interaction
{
	// point of intersection and slight offset on intersection point
	// due to precision concern
	Point3f p;
	Vector3f pError;
	// negative direction ray 
	Vector3f wo;
	// normal
	Normal3f n;
	Float time;

	BSDF* bsdf = nullptr;
	BSSRDF* bssrdf = nullptr;

	bool surfaceInteraction = false;
	
	//Interactions also need to record the scattering media at their point (if any); this is handled
	// by an instance of the MediumInterface class,
	MediumInterface mediumInterface;

	Interaction() : time(0) {}

	Interaction(const Point3f& p,
		const Normal3f& n,
		const Vector3f& pError,
		const Vector3f& wo,
		Float time,
		const MediumInterface mediumInterface);

	bool IsSurfaceInteraction() const;

};

class SurfaceInteraction : public Interaction
{
public:

	// surface u, v coordinates
	Point2f uv;
	// partial derivatives of p and n on u and v
	// they lays on the tangent plane cuz derivatives!
	Vector3f dpdu, dpdv;
	Normal3f dndu, dndv;

	const Shape* shape = nullptr;
	const Primitive* primitive = nullptr;

	// shading struct are for secondary normal values, used for stuffs like bump mappings 
	// or per-vertex normals... etc
	struct {
		Normal3f n;
		Vector3f dpdu, dpdv;
		Normal3f dndu, dndv;
	} shading;

public:

	SurfaceInteraction() {}

	SurfaceInteraction(const Point3f& p,
		const Vector3f& pError,
		const Point2f& uv,
		const Vector3f& wo,
		const Vector3f& dpdu,
		const Vector3f& dpdv,
		const Normal3f& dndu,
		const Normal3f& dndv,
		Float time,
		const Shape* shape);

	~SurfaceInteraction() {}

	// update shading struct
	void SetShadingGeometry(const Vector3f& dpdus, const Vector3f& dpdvs,
		const Normal3f& dndus, const Normal3f& dndvs, bool orientationIsAuthoritative);

	void ComputeScatteringFunctions(Ray r, MemoryArena m);
	
	Spectrum Le(Vector3f wo);



};

