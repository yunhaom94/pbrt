// flags
#pragma once
#define _USE_MATH_DEFINES

// includes
// std libs
#include <math.h>
#include <iostream>
#include <vector>
#include <memory>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>


// properties
#define PBRT_CONSTEXPR constexpr


//  Forward Declarations

template <typename T> class Bounds2;
template <typename T> class Bounds3;

template <int nSpectrumSamples> class CoefficientSpectrum;
class RGBSpectrum;
class SampledSpectrum;

class Film;

class Light;

class Options;

class Primitive;

class Ray;

class Transform;
class AnimatedTransform;

class Medium;

class RayDifferential;

template <typename T> class Texture;

class Shape;
class TriangleMesh;

class SurfaceInteraction;

class Scene;

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

// if change type of spectrum, need recompile
// may change later
typedef RGBSpectrum Spectrum;
//typedef SampledSpectrum Spectrum;


// values
static PBRT_CONSTEXPR double Pi = 3.14159265358979323846;