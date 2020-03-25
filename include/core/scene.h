#include <vector>
#include <memory>
#include "light.h"
#include "primitive.h"
#include "aggregate.h"

class Scene
{
public:
	// TODO: thing
	Scene(std::shared_ptr<Primitive> aggregate, const std::vector<std::shared_ptr<Light>>& lights) : lights(lights), aggregate(aggregate);
	~Scene();
	void Intersect();

private:

};

