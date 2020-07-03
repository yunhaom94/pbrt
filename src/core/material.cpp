#include "core/material.h"

void Material::Bump(const std::shared_ptr<Texture<Float>>& d,
	SurfaceInteraction* si) const
{
	SurfaceInteraction siEval = *si;

	Float du = 0.5 * (std::abs(si->dudx) + std::abs(si->dudy));
	if (du == 0) 
		du = 0.01;
	siEval.p = si->p + du * si->shading.dpdu;
	siEval.uv = si->uv + Vector2f(du, 0);
	siEval.n = ((Normal3f)si->shading.dpdu.cross(si->shading.dpdv) +
		du * si->dndu).normalized();

	Float uDisplace = d->Evaluate(siEval);

	Float dv = 0.5 * (std::abs(si->dvdx) + std::abs(si->dvdy));
	if (dv == 0) 
		dv = 0.01;
	siEval.p = si->p + dv * si->shading.dpdv;
	siEval.uv = si->uv + Vector2f(dv, 0);
	siEval.n = ((Normal3f)si->shading.dpdu.cross(si->shading.dpdv) +
		du * si->dndv).normalized();

	Float vDisplace = d->Evaluate(siEval);
	
	Float displace = d->Evaluate(*si);
	Vector3f dpdu = si->shading.dpdu +
		(uDisplace - displace) / du * Vector3f(si->shading.n) +
		displace * Vector3f(si->shading.dndu);
	Vector3f dpdv = si->shading.dpdv +
		(vDisplace - displace) / dv * Vector3f(si->shading.n) +
		displace * Vector3f(si->shading.dndv);
	si->SetShadingGeometry(dpdu, dpdv, si->shading.dndu, si->shading.dndv,
		false);
}
