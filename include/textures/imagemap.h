#pragma once

#include "core/texture.h"
#include "utlis/utlis.h"
#include <map>

// TODO: p618

/*
enum class ImageWrap { Repeat, Black, Clamp };

struct ResampleWeight {
    int firstTexel;
    Float weight[4];
};

// stores the texels in memory and handles the details of
// reconstruction and filtering to reduce aliasing.
template <typename T>
class MIPMap
{
private:
    const bool doTrilinear;
    const Float maxAnisotropy;
    const ImageWrap wrapMode;
    Point2i resolution;
    std::vector<std::unique_ptr<BlockedArray<T>>> pyramid;

public:
    MIPMap(const Point2i& res, const T* img, bool doTrilinear, Float maxAnisotropy, ImageWrap wrapMode);

    const T& Texel(int level, int s, int t) const;

    ~MIPMap() {}

    int Width() const { return resolution[0]; }
    int Height() const { return resolution[1]; }
    int Levels() const { return pyramid.size(); }



private:
    std::unique_ptr<ResampleWeight[]> resampleWeights(int oldRes,
        int newRes);
};


struct TexInfo
{
    std::string filename;
    bool doTrilinear;
    Float maxAniso;
    ImageWrap wrapMode;
    Float scale;
    bool gamma;

    TexInfo(const std::string& f, bool dt, Float ma, ImageWrap wm, Float sc,
        bool gamma)
        : filename(f),
        doTrilinear(dt),
        maxAniso(ma),
        wrapMode(wm),
        scale(sc),
        gamma(gamma) {}

    bool operator<(const TexInfo& t2) const
    {
        if (filename != t2.filename) return filename < t2.filename;
        if (doTrilinear != t2.doTrilinear) return doTrilinear < t2.doTrilinear;
        if (maxAniso != t2.maxAniso) return maxAniso < t2.maxAniso;
        if (scale != t2.scale) return scale < t2.scale;
        if (gamma != t2.gamma) return !gamma;
        return wrapMode < t2.wrapMode;
    }
};

template <typename Tmemory, typename Treturn>
class ImageTexture : public Texture<Treturn> 
{
private:
	std::unique_ptr<TextureMapping2D> mapping;
	MIPMap<Tmemory>* mipmap;
	static std::map<TexInfo, std::unique_ptr<MIPMap<Tmemory>>> textures;

public:
	ImageTexture(std::unique_ptr<TextureMapping2D> mapping, const std::string& filename,
		bool doTrilinear, Float maxAniso, ImageWrap wrapMode, Float scale, bool gamma);

    MIPMap<Tmemory>* GetTexture(const std::string& filename, bool doTrilinear, Float maxAniso, ImageWrap wrap, Float scale, bool gamma);

    static void ClearCache() 
    {
        textures.erase(textures.begin(), textures.end());
    }

    Treturn Evaluate(const SurfaceInteraction& si) const;



private:
    static void convertIn(const RGBSpectrum& from, RGBSpectrum* to,
        Float scale, bool gamma)
    {
        for (int i = 0; i < RGBSpectrum::nSamples; ++i)
            (*to)[i] = scale * (gamma ? InverseGammaCorrect(from[i])
                : from[i]);
    }

    static void convertIn(const RGBSpectrum& from, Float* to,
        Float scale, bool gamma)
    {
        *to = scale * (gamma ? InverseGammaCorrect(from.y())
            : from.y());
    }

    static void convertOut(const RGBSpectrum& from, Spectrum* to) 
    {
        Float rgb[3];
        from.ToRGB(rgb);
        *to = Spectrum::FromRGB(rgb);
    }

    static void convertOut(Float from, Float* to)
    {
        *to = from;
    }


};



Float Lanczos(Float, Float tau = 2);
*/