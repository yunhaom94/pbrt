#pragma once


#include "core/pbrt.h"

// TODO: p1074
class MemoryArena
{
public:
	MemoryArena() {}
	MemoryArena(int i) {}
	~MemoryArena() {}

	void Reset() {}

	template <class T>
	T* Alloc() { return nullptr; }

private:

};

