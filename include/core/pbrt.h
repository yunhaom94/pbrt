// flags
#pragma once
#define _USE_MATH_DEFINES

// includes
#include <math.h>


//  Forward Declarations and definations

class SurfaceInteraction;
class Ray;
template <typename T> class Bounds2;
template <typename T> class Bounds3;
class Transform;


// helper functions
template <typename T, typename U, typename V>
inline T Clamp(T val, U low, V high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}

// definitions
typedef Bounds2<double> Bounds2d;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<double> Bounds3d;
typedef Bounds3<int> Bounds3i;