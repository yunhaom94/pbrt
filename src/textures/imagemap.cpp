#include "textures/imagemap.h"

/*
Float Lanczos(Float x, Float tau) {
	x = std::abs(x);
	if (x < 1e-5f) return 1;
	if (x > 1.f) return 0;
	x *= Pi;
	Float s = std::sin(x * tau) / (x * tau);
	Float lanczos = std::sin(x) / x;
	return s * lanczos;
}

template <typename T>
MIPMap<T>::MIPMap(const Point2i& res, 
	const T* img,
	bool doTrilinear,
	Float maxAnisotropy, 
	ImageWrap wrapMode)
	: doTrilinear(doTrilinear), 
	maxAnisotropy(maxAnisotropy),
	wrapMode(wrapMode), 
	resolution(res)
{
	std::unique_ptr<T[]> resampledImage = nullptr;
	if (!IsPowerOf2(resolution[0]) || !IsPowerOf2(resolution[1])) 
	{
		Point2i resPow2(RoundUpPow2(resolution[0]), RoundUpPow2(resolution[1]));

	        // Resample image in $s$ direction
        std::unique_ptr<ResampleWeight[]> sWeights =
            resampleWeights(resolution[0], resPow2[0]);
        resampledImage.reset(new T[resPow2[0] * resPow2[1]]);

        // Apply _sWeights_ to zoom in $s$ direction
        ParallelFor([&](int t) {
            for (int s = 0; s < resPow2[0]; ++s) {
                // Compute texel $(s,t)$ in $s$-zoomed image
                resampledImage[t * resPow2[0] + s] = 0.f;
                for (int j = 0; j < 4; ++j) {
                    int origS = sWeights[s].firstTexel + j;
                    if (wrapMode == ImageWrap::Repeat)
                        origS = Mod(origS, resolution[0]);
                    else if (wrapMode == ImageWrap::Clamp)
                        origS = Clamp(origS, 0, resolution[0] - 1);
                    if (origS >= 0 && origS < (int)resolution[0])
                        resampledImage[t * resPow2[0] + s] +=
                            sWeights[s].weight[j] *
                            img[t * resolution[0] + origS];
                }
            }
        }, resolution[1], 16);

        // Resample image in $t$ direction
        std::unique_ptr<ResampleWeight[]> tWeights =
            resampleWeights(resolution[1], resPow2[1]);
        std::vector<T *> resampleBufs;
        int nThreads = MaxThreadIndex();
        for (int i = 0; i < nThreads; ++i)
            resampleBufs.push_back(new T[resPow2[1]]);
        ParallelFor([&](int s) {
            T *workData = resampleBufs[ThreadIndex];
            for (int t = 0; t < resPow2[1]; ++t) {
                workData[t] = 0.f;
                for (int j = 0; j < 4; ++j) {
                    int offset = tWeights[t].firstTexel + j;
                    if (wrapMode == ImageWrap::Repeat)
                        offset = Mod(offset, resolution[1]);
                    else if (wrapMode == ImageWrap::Clamp)
                        offset = Clamp(offset, 0, (int)resolution[1] - 1);
                    if (offset >= 0 && offset < (int)resolution[1])
                        workData[t] += tWeights[t].weight[j] *
                                       resampledImage[offset * resPow2[0] + s];
                }
            }
            for (int t = 0; t < resPow2[1]; ++t)
                resampledImage[t * resPow2[0] + s] = clamp(workData[t]);
        }, resPow2[0], 32);
        for (auto ptr : resampleBufs) delete[] ptr;
        resolution = resPow2;

		resolution = resPow2;
	}

	int nLevels = 1 + Log2Int(std::max(resolution[0], resolution[1]));
	pyramid.resize(nLevels);
	pyramid[0].reset(new BlockedArray<T>(resolution[0], resolution[1],
		resampledImage ? resampledImage.get() : img));

		for (int i = 1; i < nLevels; ++i) {
			int sRes = std::max(1, pyramid[i - 1]->uSize() / 2);
			int tRes = std::max(1, pyramid[i - 1]->vSize() / 2);
			pyramid[i].reset(new BlockedArray<T>(sRes, tRes));
			ParallelFor(
				[&](int t) {
				for (int s = 0; s < sRes; ++s)
					(*pyramid[i])(s, t) = .25f *
					(Texel(i - 1, 2 * s, 2 * t) + Texel(i - 1, 2 * s + 1, 2 * t) +
						Texel(i - 1, 2 * s, 2 * t + 1) + Texel(i - 1, 2 * s + 1, 2 * t + 1));
			}, tRes, 16);
		}
		Initialize EWA filter weights if needed 639
}

template <typename T>
const T& MIPMap<T>::Texel(int level, int s, int t) const {
	const BlockedArray<T>& l = *pyramid[level];
	switch (wrapMode) {
	case ImageWrap::Repeat:
		s = Mod(s, l.uSize());
		t = Mod(t, l.vSize());
		break;
	case ImageWrap::Clamp:
		s = Clamp(s, 0, l.uSize() - 1);
		t = Clamp(t, 0, l.vSize() - 1);
		break;
	case ImageWrap::Black: {
		static const T black = 0.f;
		if (s < 0 || s >= (int)l.uSize() ||
			t < 0 || t >= (int)l.vSize())
			return black;
		break;
	}
	}
		return l(s, t);
}

template<typename T>
std::unique_ptr<ResampleWeight[]> MIPMap<T>::resampleWeights(int oldRes, int newRes)
{
	Assert(newRes >= oldRes);
	std::unique_ptr<ResampleWeight[]> wt(new ResampleWeight[newRes]);
	Float filterwidth = 2.f;
	for (int i = 0; i < newRes; ++i) {
		Float center = (i + .5f) * oldRes / newRes;
		wt[i].firstTexel = std::floor((center - filterwidth) + 0.5f);
		for (int j = 0; j < 4; ++j) {
			Float pos = wt[i].firstTexel + j + .5f;
			wt[i].weight[j] = Lanczos((pos - center) / filterwidth);
		}
	}

	Float invSumWts = 1 / (wt[i].weight[0] + wt[i].weight[1] +
		wt[i].weight[2] + wt[i].weight[3]);
	for (int j = 0; j < 4; ++j)
		wt[i].weight[j] *= invSumWts;

	return wt;
}
template <typename Tmemory, typename Treturn>
ImageTexture<Tmemory, Treturn>::ImageTexture(
	std::unique_ptr<TextureMapping2D> mapping,
	const std::string& filename,
	bool doTrilinear,
	Float maxAniso,
	ImageWrap wrapMode,
	Float scale,
	bool gamma)
	: mapping(std::move(mapping))
{
	mipmap = GetTexture(filename, doTrilinear, maxAniso,
		wrapMode, scale, gamma);
}

template <typename Tmemory, typename Treturn> 
MIPMap<Tmemory>* ImageTexture<Tmemory, Treturn>::GetTexture(const std::string& filename,
	bool doTrilinear, Float maxAniso, ImageWrap wrap, Float scale,
	bool gamma)
{
	TexInfo texInfo(filename, doTrilinear, maxAniso, wrap, scale, gamma);
	if (textures.find(texInfo) != textures.end())
		return textures[texInfo].get();

	Point2i resolution;
	std::unique_ptr<RGBSpectrum[]> texels = ReadImage(filename, &resolution);
	MIPMap<Tmemory>* mipmap = nullptr;
	if (texels) {
		std::unique_ptr<Tmemory[]> convertedTexels(new Tmemory[resolution.x *
			resolution.y]);
		for (int i = 0; i < resolution.x * resolution.y; ++i)
			convertIn(texels[i], &convertedTexels[i], scale, gamma);
		mipmap = new MIPMap<Tmemory>(resolution, convertedTexels.get(),
			doTrilinear, maxAniso, wrap);
	}
	else {
		Tmemory oneVal = scale;
		mipmap = new MIPMap<Tmemory>(Point2i(1, 1), &oneVal);
	}
	textures[texInfo].reset(mipmap);

	return mipmap;
}

template<typename Tmemory, typename Treturn>
void ImageTexture<Tmemory, Treturn>::ClearCache()
{
	textures.erase(textures.begin(), textures.end());
}

template<typename Tmemory, typename Treturn>
Treturn ImageTexture<Tmemory, Treturn>::Evaluate(const SurfaceInteraction& si) const
{
	Vector2f dstdx, dstdy;
	Point2f st = mapping->Map(si, &dstdx, &dstdy);
	Tmemory mem = mipmap->Lookup(st, dstdx, dstdy);
	Treturn ret;
	convertOut(mem, &ret);
	return ret;
}
*/