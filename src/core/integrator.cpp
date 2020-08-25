#include "core/integrator.h"
#include "core/spectrum.h"
#include "core/sampler.h"
#include "core/scene.h"
#include "core/reflection.h"
#include "core/sampling.h"
#include "core/film.h"

Spectrum UniformSampleAllLights(const Interaction& it, const Scene& scene,
    MemoryArena& arena, Sampler& sampler,
    const std::vector<int>& nLightSamples, bool handleMedia)
{
    Spectrum L(0.f);
    for (size_t j = 0; j < scene.lights.size(); ++j) {
        const std::shared_ptr<Light>& light = scene.lights[j];
        int nSamples = nLightSamples[j];
        const Point2f* uLightArray = sampler.Get2DArray(nSamples);
        const Point2f* uScatteringArray = sampler.Get2DArray(nSamples);
        if (!uLightArray || !uScatteringArray) {
            Point2f uLight = sampler.Get2D();
            Point2f uScattering = sampler.Get2D();
            L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
                arena, handleMedia);
        }
        else {
            Spectrum Ld(0.f);
            for (int k = 0; k < nSamples; ++k)
                Ld += EstimateDirect(it, uScatteringArray[k], *light, uLightArray[k],
                    scene, sampler, arena, handleMedia);
            L += Ld / nSamples;
        }
    }
    return L;
}

Spectrum UniformSampleOneLight(const Interaction& it,
    const Scene& scene, MemoryArena& arena, Sampler& sampler,
    bool handleMedia)
{
    int nLights = int(scene.lights.size());
    if (nLights == 0) 
        return Spectrum(0.f);

    int lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
    const std::shared_ptr<Light>& light = scene.lights[lightNum];
    Point2f uLight = sampler.Get2D();
    Point2f uScattering = sampler.Get2D();

    return (Float)nLights *
        EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
            arena, handleMedia);
}


Spectrum EstimateDirect(const Interaction& it,
    const Point2f& uScattering, const Light& light,
    const Point2f& uLight, const Scene& scene, Sampler& sampler,
    MemoryArena& arena, bool handleMedia, bool specular)
{
    BxDFType bsdfFlags = specular ? BSDF_ALL :
        BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
    Spectrum Ld(0.f);
    
    //Sample light source with multiple importance sampling 
    Vector3f wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;
    Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);

    if (lightPdf > 0 && !Li.IsBlack()) 
    {
        Spectrum f;
        if (it.IsSurfaceInteraction())
        {
            if (handleMedia)
                Li *= visibility.Tr(scene, sampler);
            else if (!visibility.Unoccluded(scene))
                Li = Spectrum(0);
        }
        else
        {
            // TODO: Evaluate phase function for light sampling strategy 900
        }
            
        if (!f.IsBlack())
        {
            if (handleMedia)
                Li *= visibility.Tr(scene, sampler);
            else if (!visibility.Unoccluded(scene))
                Li = Spectrum(0.f);
           
            if (!Li.IsBlack()) {
                if (IsDeltaLight(light.flags))
                    Ld += f * Li / lightPdf;
                else {
                    Float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
                    Ld += f * Li * weight / lightPdf;
                }
            }
        }
    }
    
    //Sample BSDF with multiple importance sampling 860
    if (!IsDeltaLight(light.flags)) 
    {
        Spectrum f;
        bool sampledSpecular = false;
        if (it.IsSurfaceInteraction())
        {
            BxDFType sampledType;
            const SurfaceInteraction& isect = (const SurfaceInteraction&)it;
            f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
                bsdfFlags, &sampledType);
            f *= std::abs(wi.dot(isect.shading.n));
            sampledSpecular = sampledType & BSDF_SPECULAR;
        }
        else 
        {
            // TODO: Sample scattered direction for medium interactions 900
        }
        if (!f.IsBlack() && scatteringPdf > 0)
        {
            Float weight = 1;
            if (!sampledSpecular) {
                lightPdf = light.Pdf_Li(it, wi);
                if (lightPdf == 0)
                    return Ld;
                weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
            }
            
            SurfaceInteraction lightIsect;
            Ray ray = it.SpawnRay(wi);
            Spectrum Tr(1.f);
            bool foundSurfaceInteraction = handleMedia ?
                scene.IntersectTr(ray, sampler, &lightIsect, &Tr) :
                scene.Intersect(ray, &lightIsect);

            Spectrum Li(0.f);
            if (foundSurfaceInteraction) {
                if (lightIsect.primitive->GetAreaLight() == &light)
                    Li = lightIsect.Le(-wi);
            }
            else
                Li = light.Le(ray);
            if (!Li.IsBlack())
                Ld += f * Li * Tr * weight / scatteringPdf;
        }
    }
       
    return Ld;
}


void SamplerIntegrator::Render(const Scene& scene)
{
    Preprocess(scene, *sampler);

    Bounds2i sampleBounds = camera->film->GetSampleBounds();
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x() + tileSize - 1) / tileSize,
        (sampleExtent.y() + tileSize - 1) / tileSize);

    //TODO: ParallelFor2D Render image tiles in parallel p26 p1088
    for (Point2i tile = Point2i(0, 0); 1 < 1;)
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

        // loop through each pixel in the tile
        for (Point2i pixel : tileBounds) 
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
    camera->film->WriteImage(1);

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
