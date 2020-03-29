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

	// reverse the direction of normals if needed
	// so that they face the same direction 
	if (orientationIsAuthoritative)
		n = n.dot(shading.n) < 0 ? -n : n;
	else
		shading.n = shading.n.dot(n) < 0 ? -shading.n : shading.n; 

	shading.dpdu = dpdus;
	shading.dpdv = dpdvs;
	shading.dndu = dndus;
	shading.dndv = dndvs;

}
