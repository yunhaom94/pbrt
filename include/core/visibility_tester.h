#pragma once

#include "pbrt.h"

// TODO: p717
class VisibilityTester
{
public:
	VisibilityTester() {}
	~VisibilityTester() {}

	virtual bool Unoccluded(Scene s) {}

private:

};

