#pragma once


#include "core/pbrt.h"
#include <list>
#include <cstddef>

// TODO: p1074
class MemoryArena
{
public:
	MemoryArena() {}
	MemoryArena(int i) {}
	~MemoryArena() {}

	void Reset() {}


	template <typename T>
	T* Alloc(size_t size = 1) { return nullptr; }

	void *Alloc(size_t size = 1) { return nullptr; }

private:

};

// Memory Declarations
#define ARENA_ALLOC(arena, Type) new ((arena).Alloc(sizeof(Type))) Type

void* AllocAligned(size_t size);

template <typename T>
T* AllocAligned(size_t count)
{
	return (T*)AllocAligned(count * sizeof(T));
}

inline void* AllocAligned(size_t size)
{
	return _aligned_malloc(size, PBRT_L1_CACHE_LINE_SIZE);
}

inline void FreeAligned(void* ptr)
{
	if (!ptr) return;

	_aligned_free(ptr);
}