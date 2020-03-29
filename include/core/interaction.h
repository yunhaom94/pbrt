#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "core/medium_interface.h"
#include "shapes/shape.h"


struct Interaction
{
	// point of intersection and slight offset on intersection point
	// due to precision concern
	Eigen::Vector3d p, pError;
	// negative ray
	Eigen::Vector3d wo;
	// normal
	Eigen::Vector3d n;
	double time;
	
	bool surfaceInteraction = false;
	
	//Interactions also need to record the scattering media at their point (if any); this is handled
	// by an instance of the MediumInterface class,
	MediumInterface mediumInterface;

	Interaction() {}

	Interaction(const  Eigen::Vector3d &p, 
		const Eigen::Vector3d &n, 
		const Eigen::Vector3d &pError,
		const Eigen::Vector3d &wo, 
		double time,
		const MediumInterface &mediumInterface)
		: p(p), pError(pError), time(time), wo(wo), n(n),
		mediumInterface(mediumInterface) { }


	// TODO: p115
	bool IsSurfaceInteraction() const;

};

class SurfaceInteraction : public Interaction
{
public:

	// surface coordinates
	Eigen::Vector2d uv;
	// partial derivatives of p and n on u and v
	// they lays on the tangent plane cuz derivatives!
	Eigen::Vector3d dpdu, dpdv, dndu, dndv;

	const Shape* shape = nullptr;

	// shading struct are for secondary normal values, used for stuffs like bump mappings 
	// or per-vertex nromals... etc
	struct {
		Eigen::Vector3d n, dpdu, dpdv, dndu, dndv;
	} shading;

public:
	SurfaceInteraction() {}

	SurfaceInteraction(const Eigen::Vector3d& p,
		const Eigen::Vector3d& pError,
		const Eigen::Vector2d& uv,
		const Eigen::Vector3d& wo,
		const Eigen::Vector3d& dpdu,
		const Eigen::Vector3d& dpdv,
		const Eigen::Vector3d& dndu,
		const Eigen::Vector3d& dndv,
		double time,
		const Shape* shape)
		: Interaction(p, dpdu.cross(dpdv).normalized(), pError, wo, time, MediumInterface()),
		uv(uv), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv), shape(shape)
	{
		surfaceInteraction = true;
		shading.n = n;
		shading.dpdu = dpdu;
		shading.dpdv = dpdv;
		shading.dndu = dndu;
		shading.dndv = dndv;

		// sometime normals needs to be reversed (specified in the input)
		if (shape && (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
			n *= -1;
			shading.n *= -1;
		}
	}

	// update shading struct
	void SetShadingGeometry(const Eigen::Vector3d& dpdus, const Eigen::Vector3d& dpdvs,
		const Eigen::Vector3d& dndus, const Eigen::Vector3d& dndvs, bool orientationIsAuthoritative);

	~SurfaceInteraction() {}

private:

};

