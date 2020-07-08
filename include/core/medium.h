#pragma once
class Medium
{
public:
	Medium();
	~Medium();

	virtual Spectrum Tr(const Ray& ray, Sampler& sampler) const = 0;

private:

};

