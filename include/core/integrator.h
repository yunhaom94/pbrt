#include "scene.h"

class Integrator
{
public:
	Integrator();
	~Integrator();
	virtual void Render(const Scene& scene)) = 0;

private:

};

Integrator::Integrator()
{
}

Integrator::~Integrator()
{
}