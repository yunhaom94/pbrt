#pragma once
#include "core/pbrt.h"
#include "core/bounding_boxes.h"
#include "core/filter.h"
#include "core/spectrum.h"
#include <mutex>
#include <fstream>

// TODO: p1068
inline void WriteImage(const std::string& name, const Float* rgb,
	const Bounds2i& outputBounds, const Point2i& totalResolution)
{
	std::ofstream outfile(name);
	int image_height = totalResolution.y();
	int image_width = totalResolution.x();
	outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";


	for (int j = 0; j < image_height; j++)
	{
		std::cerr << "\rScanlines remaining: " << image_height - j - 1 << ' ' << std::flush;
		for (int i = 0; i < image_width; i++)
		{
			int offset = j * (3 * image_width) + i * 3;
			Float color[3] = { rgb[offset], rgb[offset + 1], rgb[offset + 2] };

			outfile << (int)Clamp<Float, Float, Float>(255 * color[0], 0, 255.0) << ' '
				<< (int)Clamp<Float, Float, Float>(255 * color[1], 0, 255.0) << ' '
				<< (int)Clamp<Float, Float, Float>(255 * color[2], 0, 255.0) << '\n';
		}
	}

}

// for each pixel
struct FilmTilePixel
{
	// sum of the weighted contributions from the pixel samples
	Spectrum contribSum = 0.0;
	// sum of the filter weights is maintained.
	Float filterWeightSum = 0.0;
};

// Stores contributions for the pixels 
// in the corresponding region of the image.
// Useful in multi-threading.
class FilmTile
{
private:
	// in x * y pixel numbers
	const Bounds2i pixelBounds;
	const Vector2f filterRadius, invFilterRadius;
	const Float* filterTable;
	const int filterTableSize;
	std::vector<FilmTilePixel> pixels;

public:
	FilmTile(const Bounds2i& pixelBounds, const Vector2f& filterRadius, const Float* filterTable, int filterTableSize);

	// add a sample value from a ray
	void AddSample(const Point2f& pFilm, const Spectrum& L,
		Float sampleWeight = 1);

	FilmTilePixel& GetPixel(const Point2i& p);

	Bounds2i GetPixelBounds() const { return pixelBounds; }


};

// film is the "imagine plane" that camera shoots ray to
// Film class represent the pixel definition of a final image
class Film
{
public:
	// resolution of the image in pixels;
	const Point2i fullResolution;
	// length of the diagonal of screen(film) in meters
	const Float diagonal;
	std::unique_ptr<Filter> filter;
	// output filename
	const std::string filename;
	// a crop window bound that converted into pixels
	Bounds2i croppedPixelBounds;

private:
	struct Pixel
	{	
		// weighted sums of XYZ colors
		Float xyz[3] = { 0, 0, 0 };
		// sum of filter weight values for the sample
		// contributions
		Float filterWeightSum = 0;
		// (unweighted) sum of sample splats.
		AtomicFloat splatXYZ[3] = { 0, 0, 0};
		Float pad; //unused for padding
	};
	
	// frame buffer
	std::unique_ptr<Pixel[]> pixels;
	// # of pixels that image sample contribute to
	static constexpr int filterTableWidth = 16;
	// precomputation a table of filter values, basically a kernel (?)
	Float filterTable[filterTableWidth * filterTableWidth];
	std::mutex mutex;
	Float scale;

public:
	// cropWindow: a crop window that may specify a subset of the image to render in [0-1]
	// Diagonal: length of the diagonal of physical film in millimeters and converted to meters
	Film(const Point2i& resolution, const Bounds2f& cropWindow,
		std::unique_ptr<Filter> filt, Float diagonal,
		const std::string& filename, Float scale);

	~Film() {}

	// give the sampler a bound slightly larger then the film bounds
	// so that edges don't loss info on filtering
	Bounds2i GetSampleBounds() const;

	// used by Realistic Camera
	Bounds2f GetPhysicalExtent() const;

	// get a region on film, used for multi-threading
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


