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
class FilmTile;

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

class Sampler;

class Camera;

class MemoryArena;

// definitions
typedef double Float;

typedef Bounds2<double> Bounds2d;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<double> Bounds3d;
typedef Bounds3<int> Bounds3i;

// if change type of spectrum, need recompile
// may change later
typedef RGBSpectrum Spectrum;
//typedef SampledSpectrum Spectrum;

// define eigen vectors to match pbrt vectors
typedef Eigen::Vector2i Vector2i;
typedef Eigen::Vector2i Point2i;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector2d Point2d;

// values
static PBRT_CONSTEXPR double Pi = 3.14159265358979323846;
