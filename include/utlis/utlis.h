#pragma once
#include "core/pbrt.h"

// linearly interpolate
inline double Lerp(double v0, double v1, double t) {
    return (1 - t) * v0 + t * v1;
}

template <typename Predicate> 
inline int FindInterval(int size,
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

template <typename T, typename U, typename V>
inline T Clamp(T val, U low, V high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}

template <typename T> 
inline Eigen::Matrix<T, 3, 1> Faceforward(const Eigen::Matrix<T, 3, 1>& n, const Eigen::Matrix<T, 3, 1>& v) {
        return (n.dot(v) < 0) ? -n : n;
    }