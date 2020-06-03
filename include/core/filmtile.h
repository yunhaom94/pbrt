#pragma once

#include "core/pbrt.h"
#include "core/spectrum.h"

// TODO: p489, move implementation and extra includes to .cpp
class FilmTile
{
public:
	FilmTile();
	~FilmTile();
	void AddSample(Point2d pf, Spectrum s, Float weight) {}



private:

};

FilmTile::FilmTile()
{
}

FilmTile::~FilmTile()
{
}