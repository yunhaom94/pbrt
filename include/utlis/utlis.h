#pragma once
#include "core/pbrt.h"

// linearly interpolate
double Lerp(double v0, double v1, double t) {
    return (1 - t) * v0 + t * v1;
}

template <typename Predicate> int FindInterval(int size,
    const Predicate& pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        }
        else
            len = half;
    }
    return Clamp(first - 1, 0, size - 2);
}