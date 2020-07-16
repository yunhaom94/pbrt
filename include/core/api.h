#pragma once

#include "core/pbrt.h"

template <typename T>
struct ParamSetItem
{
	const std::string name;
	const std::unique_ptr<T[]> values;
	const int nValues;
	mutable bool lookedUp = false;

	template <typename T>
	ParamSetItem<T>::ParamSetItem(const std::string& name, const T* v,
		int nValues)
		: name(name), values(new T[nValues]), nValues(nValues) 
	{
		std::copy(v, v + nValues, values.get());
	}
};

class ParamSet
{
private:
	std::vector<std::shared_ptr<ParamSetItem<bool>>> bools;
	std::vector<std::shared_ptr<ParamSetItem<int>>> ints;
	std::vector<std::shared_ptr<ParamSetItem<Float>>> floats;
	std::vector<std::shared_ptr<ParamSetItem<Point2f>>> point2fs;
	std::vector<std::shared_ptr<ParamSetItem<Vector2f>>> vector2fs;
	std::vector<std::shared_ptr<ParamSetItem<Point3f>>> point3fs;
	std::vector<std::shared_ptr<ParamSetItem<Vector3f>>> vector3fs;
	std::vector<std::shared_ptr<ParamSetItem<Normal3f>>> normals;
	std::vector<std::shared_ptr<ParamSetItem<Spectrum>>> spectra;
	std::vector<std::shared_ptr<ParamSetItem<std::string>>> strings;
	std::vector<std::shared_ptr<ParamSetItem<std::string>>> textures;

public:

};