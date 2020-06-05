#include "integrators/sampler_integrator.h"
#include "core/film.h"
#include "core/camera.h"
#include "core/bounding_boxes.h"
#include "core/memory.h"
#include "core/sampler.h"
#include "core/filmtile.h"
#include "core/ray_differential.h"
#include "core/spectrum.h"
#include "core/interaction.h"
#include "core/bxdf.h"

void SamplerIntegrator::Render(const Scene& scene)
{
	Preprocess(scene, *sampler);

	Bounds2i sampleBounds = camera->film->GetSampleBounds();
	Vector2i sampleExtent = sampleBounds.Diagonal();
	const int tileSize = 16;
	Point2i nTiles((sampleExtent.x() + tileSize - 1) / tileSize,
		(sampleExtent.y() + tileSize - 1) / tileSize);

	//TODO: ParallelFor2D Render image tiles in parallel p26 p1088
	for (Point2i tile = Point2i(0,0); 1 < 1;)
	{
		MemoryArena arena;

		// Get sampler instance for tile
		int seed = tile.y() * nTiles.x() + tile.x();
		std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

		// Compute sample bounds for tile
		// because tile size may be out of the actual picture
		int x0 = sampleBounds.pMin.x() + tile.x() * tileSize;
		int y0 = sampleBounds.pMin.y() + tile.y() * tileSize;

		int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x());
		int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y());
		Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

		// a small buffer of memory to store pixel values for the current tile.
		std::unique_ptr<FilmTile> filmTile = camera->film->GetFilmTile(tileBounds);

		//TODO: implement iterator for bounding boxes p76 and p30
		// loop through each pixel in the tile
		//for (Point2i pixel : tileBounds) 
		for (Point2i pixel(0,0);;)
		{
			tileSampler->StartPixel(pixel);
			do {
				// get position on film that camera should shot rays at
				CameraSample cameraSample = tileSampler->GetCameraSample(pixel);

				// ray is actually the main ray and rays generated 1-pixel away on both x, y
				// on image plane, this is for better texture quality
				RayDifferential ray;
				
				// sometime each ray may have different weight
				Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);

				// scales the differential rays to account for the actual spacing 
				// between samples on the film plane for the case where multiple samples are taken per pixel.
				ray.ScaleDifferentials(1 / std::sqrt(tileSampler->samplesPerPixel));

				Spectrum L(0.0);
				if (rayWeight > 0)
					// calculate light radiance at the start the ray;
					// getting the actual color spectrum gathered from the ray
					L = Li(ray, scene, *tileSampler, arena);
					//TODO: Issue warning if unexpected radiance value is returned

				// add result to the image (plane)
				filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

				// memory thing
				arena.Reset();

			} while (tileSampler->StartNextSample());

			//transfer ownership of the unique_ptr to MergeFilmTile().
			camera->film->MergeFilmTile(std::move(filmTile));
		}
	}

	// output
	camera->film->WriteImage();

}

Spectrum SamplerIntegrator::SpecularReflect(const RayDifferential& ray, 
	const SurfaceInteraction& isect, 
	const Scene& scene, 
	Sampler& sampler, 
	MemoryArena& arena, 
	int depth) const
{
	// == Vector3f wo = insect wo; Vector3f wi; ????
	Vector3f wo = isect.wo, wi;
	Float pdf;

	// set type to relection for sampler
	BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
	// Given wo, compute specular reflection direction wi and BSDF value
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

	const Normal3f& ns = isect.shading.n;
	if (pdf > 0 && !f.IsBlack() && std::abs(wi.dot(ns)) != 0) 
	{
			//TODO: Compute ray differential rd for specular reflection p607
		Ray rd;
		return f * Li(rd, scene, sampler, arena, depth + 1) * std::abs(wi.dot(ns)) * (1 / pdf);
	}
	else
		return Spectrum(0.f);

	return Spectrum();
}

Spectrum SamplerIntegrator::SpecularTransmit(const RayDifferential& ray,
	const SurfaceInteraction& isect, 
	const Scene& scene, 
	Sampler& sampler,
	MemoryArena& arena, 
	int depth) const
{
	// == Vector3f wo = insect wo; Vector3f wi; ????
	Vector3f wo = isect.wo, wi;
	Float pdf;

	// set type to relection for sampler
	BxDFType type = BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR);
	// Given wo, compute specular reflection direction wi and BSDF value
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

	const Normal3f& ns = isect.shading.n;
	if (pdf > 0 && !f.IsBlack() && std::abs(wi.dot(ns)) != 0)
	{
		//TODO: Compute ray differential rd for specular reflection p607
		Ray rd;
		return f * Li(rd, scene, sampler, arena, depth + 1) * std::abs(wi.dot(ns)) * (1 / pdf);
	}
	else
		return Spectrum(0.f);

	return Spectrum();
}
