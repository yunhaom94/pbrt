#include "core/interaction.h"


bool Interaction::IsSurfaceInteraction() const
{
	return surfaceInteraction;
}

void SurfaceInteraction::SetShadingGeometry(const Eigen::Vector3d& dpdus, const Eigen::Vector3d& dpdvs, const Eigen::Vector3d& dndus, const Eigen::Vector3d& dndvs, bool orientationIsAuthoritative)
{
	shading.n = dpdus.cross(dpdvs).normalized();
	if (shape && (shape->reverseOrientation ^
		shape->transformSwapsHandedness))
		shading.n = -shading.n;

	// TODO: p119, p72
	//if (orientationIsAuthoritative)
		//n = Faceforward(n, shading.n);
	//else
		//shading.n = Faceforward(shading.n, n);
	shading.dpdu = dpdus;
	shading.dpdv = dpdvs;
	shading.dndu = dndus;
	shading.dndv = dndvs;

}
