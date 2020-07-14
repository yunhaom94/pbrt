#pragma once


#include "core/pbrt.h"
#include <list>
#include <cstddef>

// Memory Declarations
#define ALLOCA(TYPE, COUNT) (TYPE *)alloca((COUNT) * sizeof(TYPE))
#define ARENA_ALLOC(arena, Type) new ((arena).Alloc(sizeof(Type))) Type

void* AllocAligned(size_t size);

template <typename T>
T* AllocAligned(size_t count)
{
	return (T*)AllocAligned(count * sizeof(T));
}

void FreeAligned(void* ptr);

class MemoryArena
{
private:
	const size_t blockSize; // in # of bytes
	size_t currentBlockPos = 0, currentAllocSize = 0;
	uint8_t* currentBlock = nullptr;
	std::list<std::pair<size_t, uint8_t*>> usedBlocks, availableBlocks;

public:
	MemoryArena(size_t blockSize = 262144) : blockSize(blockSize) { }
	~MemoryArena() {}

	void* Alloc(size_t nBytes);

	template<typename T> 
	T* Alloc(size_t n = 1, bool runConstructor = true)
	{
		T* ret = (T*)Alloc(n * sizeof(T));
		if (runConstructor)
			for (size_t i = 0; i < n; ++i)
				new (&ret[i]) T();
		return ret;
	}

	void Reset() {}
};

