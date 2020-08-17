#include "core/integrator.h"
#include "core/spectrum.h"
#include "core/sampler.h"
#include "core/scene.h"
#include "core/reflection.h"
#include "core/sampling.h"

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