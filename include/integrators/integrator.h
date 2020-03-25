#pragma once
class Integrator
{
public:
	Integrator();
	~Integrator();
	virtual void Render() = 0;

private:

};

Integrator::Integrator()
{
}

Integrator::~Integrator()
{
}