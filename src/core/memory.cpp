#include "core/memory.h"


void* AllocAligned(size_t size)
{
#if defined(PBRT_IS_WINDOWS)
	return _aligned_malloc(size, PBRT_L1_CACHE_LINE_SIZE);
#elif defined (PBRT_IS_OPENBSD) || defined(PBRT_IS_OSX)
	void* ptr;
	if (posix_memalign(&ptr, PBRT_L1_CACHE_LINE_SIZE, size) != 0)
		ptr = nullptr;
	return ptr;
#else
	return memalign(PBRT_L1_CACHE_LINE_SIZE, size);
#endif
}



void FreeAligned(void* ptr)
{
	if (!ptr)
		return;

	_aligned_free(ptr);
}

void* MemoryArena::Alloc(size_t nBytes)
{
	nBytes = ((nBytes + 15) & (~15));

	if (currentBlockPos + nBytes > currentAllocSize)
	{
		if (currentBlock) 
		{
			usedBlocks.push_back(std::make_pair(currentAllocSize, currentBlock));
			currentBlock = nullptr;
		}
		for (auto iter = availableBlocks.begin(); iter != availableBlocks.end();
			++iter) 
		{
			if (iter->first >= nBytes) 
			{
				currentAllocSize = iter->first;
				currentBlock = iter->second;
				availableBlocks.erase(iter);
				break;
			}

			if (!currentBlock) 
			{
				currentAllocSize = std::max(nBytes, blockSize);
				currentBlock = AllocAligned<uint8_t>(currentAllocSize);
			}
			currentBlockPos = 0;

		}
	}
	void* ret = currentBlock + currentBlockPos;
	currentBlockPos += nBytes;
	return ret;
}

