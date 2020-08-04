#pragma once
#include "core/pbrt.h"
#include "core/bounding_boxes.h"
#include "core/filter.h"
#include "core/spectrum.h"
#include <mutex>

// TODO: p1068
inline void WriteImage(const std::string& name, const Float* rgb,
	const Bounds2i& outputBounds, const Point2i& totalResolution)
{

}

// film is the "imagine plane" that camera shoots ray to
class Film
{
public:
	// resolution of the image in pixels;
	const Point2i fullResolution;
	const Float diagonal;
	std::unique_ptr<Filter> filter;
	// TODO: maybe write to file with another class
	const std::string filename;
	Bounds2i croppedPixelBounds;

private:
	struct Pixel
	{
		Float xyz[3] = { 0, 0, 0 };
		Float filterWeightSum = 0;
		AtomicFloat splatXYZ[3];
		//Float pad;
	};
	std::unique_ptr<Pixel[]> pixels;
	static constexpr int filterTableWidth = 16;
	Float filterTable[filterTableWidth * filterTableWidth];
	std::mutex mutex;
	Float scale;

public:
	// cropWindow: a crop window that may specify a subset of the image to render in [0-1]
	Film(const Point2i& resolution, const Bounds2f& cropWindow,
		std::unique_ptr<Filter> filt, Float diagonal,
		const std::string& filename, Float scale);

	~Film() {}

	Bounds2i GetSampleBounds() const;

	Bounds2f GetPhysicalExtent() const;

	std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i& sampleBounds);

	void MergeFilmTile(std::unique_ptr<FilmTile> tile);

	void SetImage(const Spectrum* img) const;

	void AddSplat(const Point2f& p, const Spectrum& v);

	void WriteImage(Float splatScale);

private:
	Pixel& GetPixel(const Point2i& p)
	{
		int width = croppedPixelBounds.pMax.x() - croppedPixelBounds.pMin.x();
		int offset = (p.x() - croppedPixelBounds.pMin.x()) +
			(p.y() - croppedPixelBounds.pMin.y()) * width;
		return pixels[offset];
	}

};


struct FilmTilePixel
{
	Spectrum contribSum = 0;
	Float filterWeightSum = 0;

};

class FilmTile 
{
private:
	const Bounds2i pixelBounds;
	const Vector2f filterRadius, invFilterRadius;
	const Float* filterTable;
	const int filterTableSize;
	std::vector<FilmTilePixel> pixels;

public:
	FilmTile(const Bounds2i& pixelBounds, const Vector2f& filterRadius, const Float* filterTable, int filterTableSize);

	void AddSample(const Point2f& pFilm, const Spectrum& L,
		Float sampleWeight = 1);

	FilmTilePixel& GetPixel(const Point2i& p);

	Bounds2i GetPixelBounds() const { return pixelBounds; }


};