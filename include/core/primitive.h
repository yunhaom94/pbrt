#pragma once
class Primitive
{
public:
	Primitive();
	~Primitive();
	void virtual IntersectP() = 0;

private:

};

Primitive::Primitive()
{
}

Primitive::~Primitive()
{
}