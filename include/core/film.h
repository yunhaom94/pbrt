#pragma once
#include "core/pbrt.h"
#include "core/bounding_boxes.h"
#include "core/filmtile.h"

// TODO: p484, move implementation and extra includes to .cpp
// film is the "imagine plane" that camera shoots ray to
class Film
{
public:
	Film() {}
	~Film() {}

	Bounds2i GetSampleBounds()
	{

		return *(new Bounds2i);
	}

	std::unique_ptr<FilmTile> GetFilmTile(Bounds2i tileBounds)
	{
		return nullptr;
	}

	void MergeFilmTile(std::unique_ptr<FilmTile> t) {}

	void WriteImage() {}

private:

};

