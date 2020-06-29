#include "core/film.h"
#include "utlis/geometry.h"


Film::Film(const Point2i& resolution, const Bounds2f& cropWindow,
	std::unique_ptr<Filter> filt, Float diagonal,
	const std::string& filename, Float scale)
	: fullResolution(resolution), diagonal(diagonal * .001),
	filter(std::move(filt)), filename(filename), scale(scale)
{
	croppedPixelBounds =
		Bounds2i(Point2i(std::ceil(fullResolution.x() * cropWindow.pMin.x()),
			std::ceil(fullResolution.y() * cropWindow.pMin.y())),
			Point2i(std::ceil(fullResolution.x() * cropWindow.pMax.x()),
				std::ceil(fullResolution.y() * cropWindow.pMax.y())));

	pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);

	int offset = 0;
	for (int y = 0; y < filterTableWidth; ++y) {
		for (int x = 0; x < filterTableWidth; ++x, ++offset) {
			Point2f p;
			p.x() = (x + 0.5f) * filter->radius.x() / filterTableWidth;
			p.y() = (y + 0.5f) * filter->radius.y() / filterTableWidth;
			filterTable[offset] = filter->Evaluate(p);
		}
	}
}

Bounds2i Film::GetSampleBounds() const 
{
	Point2f a = Point2f(croppedPixelBounds.pMin) + Point2f(0.5, 0.5) - (Point2f)filter->radius;

	Point2f b = Point2f(croppedPixelBounds.pMax) - Point2f(0.5, 0.5) + (Point2f)filter->radius;

	Bounds2f floatBounds(Floor(a), Ceil(b));

	return (Bounds2i)floatBounds;
}

Bounds2f Film::GetPhysicalExtent() const 
{
	Float aspect = (Float)fullResolution.y() / (Float)fullResolution.x();
	Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
	Float y = aspect * x;
	return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

std::unique_ptr<FilmTile> Film::GetFilmTile(
	const Bounds2i& sampleBounds)
{
	Vector2f halfPixel = Vector2f(0.5f, 0.5f);
	Bounds2f floatBounds = (Bounds2f)sampleBounds;

	Point2i p0 = (Point2i)Ceil(Point2f(floatBounds.pMin - 
		(Point2f)halfPixel - 
		(Point2f)filter->radius));

	Point2i p1 = (Point2i)Floor(Point2f(floatBounds.pMax -
		(Point2f)halfPixel + 
		(Point2f)filter->radius)) + Point2i(1, 1);

	Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);

	return std::unique_ptr<FilmTile>(new FilmTile(tilePixelBounds,
			filter->radius, filterTable, filterTableWidth));
}



void FilmTile::AddSample(const Point2f& pFilm, const Spectrum& L, Float sampleWeight)
{
	Point2f pFilmDiscrete = pFilm - Vector2f(0.5f, 0.5f);
	Point2i p0 = (Point2i)Ceil(Point2f(pFilmDiscrete - filterRadius));
	Point2i p1 = (Point2i)Floor(Point2f(pFilmDiscrete + filterRadius)) + Point2i(1, 1);
	p0 = Max(p0, pixelBounds.pMin);
	p1 = Min(p1, pixelBounds.pMax);

	int* ifx = ALLOCA(int, p1.x() - p0.x());
	for (int x = p0.x(); x < p1.x(); ++x) {
		Float fx = std::abs((x - pFilmDiscrete.x()) *
			invFilterRadius.x() * filterTableSize);
		ifx[x - p0.x()] = std::min((int)std::floor(fx), filterTableSize - 1);
	}
	int* ify = ALLOCA(int, p1.y() - p0.y());
	for (int y = p0.y(); y < p1.y(); ++y) {
		Float fy = std::abs((y - pFilmDiscrete.y()) *
			invFilterRadius.y() * filterTableSize);
		ify[y - p0.y()] = std::min((int)std::floor(fy), filterTableSize - 1);
	}

	for (int y = p0.y(); y < p1.y(); ++y) {
		for (int x = p0.x(); x < p1.x(); ++x)
		{
			int offset = ify[y - p0.y()] * filterTableSize + ifx[x - p0.x()];
			Float filterWeight = filterTable[offset];
			FilmTilePixel& pixel = GetPixel(Point2i(x, y));
			pixel.contribSum += L * sampleWeight * filterWeight;
			pixel.filterWeightSum += filterWeight;
		}
	}
}

void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) 
{
	std::lock_guard<std::mutex> lock(mutex);
	for (Point2i pixel : tile->GetPixelBounds()) {
		const FilmTilePixel& tilePixel = tile->GetPixel(pixel);
		Pixel& mergePixel = GetPixel(pixel);
		Float xyz[3];
		tilePixel.contribSum.ToXYZ(xyz);
		for (int i = 0; i < 3; ++i)
			mergePixel.xyz[i] += xyz[i];
		mergePixel.filterWeightSum += tilePixel.filterWeightSum;
	}
}

void Film::SetImage(const Spectrum* img) const {
	int nPixels = croppedPixelBounds.Area();
	for (int i = 0; i < nPixels; ++i) {
		Pixel& p = pixels[i];
		img[i].ToXYZ(p.xyz);
		p.filterWeightSum = 1;
		p.splatXYZ[0] = p.splatXYZ[1] = p.splatXYZ[2] = 0;
	}
}

void Film::AddSplat(const Point2f& p, const Spectrum& v) {
	if (!InsideExclusive((Point2i)p, croppedPixelBounds))
		return;
	Float xyz[3];
	v.ToXYZ(xyz);
	Pixel& pixel = GetPixel((Point2i)p);
	for (int i = 0; i < 3; ++i)
		pixel.splatXYZ[i] + xyz[i];// TODO: p1087 pixel.splatXYZ[i].Add(xyz[i]);
}

void Film::WriteImage(Float splatScale)
{
	std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
	int offset = 0;
	for (Point2i p : croppedPixelBounds) 
	{
		Pixel& pixel = GetPixel(p);
		XYZToRGB(pixel.xyz, &rgb[3 * offset]);
		Float filterWeightSum = pixel.filterWeightSum;
		if (filterWeightSum != 0) {
			Float invWt = (Float)1 / filterWeightSum;
			rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
			rgb[3 * offset + 1] = std::max((Float)0, rgb[3 * offset + 1] * invWt);
			rgb[3 * offset + 2] = std::max((Float)0, rgb[3 * offset + 2] * invWt);
		}
		Float splatRGB[3];
		Float splatXYZ[3] = { pixel.splatXYZ[0], pixel.splatXYZ[1],
		pixel.splatXYZ[2] };
		XYZToRGB(splatXYZ, splatRGB);
		rgb[3 * offset] += splatScale * splatRGB[0];
		rgb[3 * offset + 1] += splatScale * splatRGB[1];
		rgb[3 * offset + 2] += splatScale * splatRGB[2];
		rgb[3 * offset] *= scale;
		rgb[3 * offset + 1] *= scale;
		rgb[3 * offset + 2] *= scale;

		++offset;
	}
	::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
}


FilmTile::FilmTile(const Bounds2i& pixelBounds,
	const Vector2f& filterRadius,
	const Float* filterTable,
	int filterTableSize)
	: pixelBounds(pixelBounds),
	filterRadius(filterRadius),
	invFilterRadius(1 / filterRadius.x(), 1 / filterRadius.y()),
	filterTable(filterTable), filterTableSize(filterTableSize)
{
	pixels = std::vector<FilmTilePixel>(std::max(0, pixelBounds.Area()));
}

FilmTilePixel& FilmTile::GetPixel(const Point2i& p) 
{
	int width = pixelBounds.pMax.x() - pixelBounds.pMin.x();
	int offset = (p.x() - pixelBounds.pMin.x()) +
		(p.y() - pixelBounds.pMin.y()) * width;
	return pixels[offset];
}

