// flags
#pragma once
#define _USE_MATH_DEFINES

// includes
// std libs
#include <math.h>
#include <iostream>
#include <vector>
#include <memory>
#include <limits>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// properties and constants
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

class VisibilityTester;

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

class Interaction;
class SurfaceInteraction;

class Scene;

class Sampler;

class Camera;

class MemoryArena;

class BSDF;

class MediumInterface;

// definitions
// TODO: may have flag for float and double
typedef double Float;

typedef Bounds2<Float> Bounds2f;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<Float> Bounds3f;
typedef Bounds3<int> Bounds3i;

// if change type of spectrum, need recompile
// may change later
typedef RGBSpectrum Spectrum;
//typedef SampledSpectrum Spectrum;

// define eigen vectors to match pbrt vectors
// TODO: may have flag for float and double
typedef Eigen::Vector2i Vector2i;
typedef Eigen::Vector2i Point2i;
typedef Eigen::Vector2d Vector2f;
typedef Eigen::Vector2d Point2f;

typedef Eigen::Vector3d Normal3f;
typedef Eigen::Vector3d Vector3f;
typedef Eigen::Vector3i Point3i;
typedef Eigen::Vector3d Point3f;

typedef Eigen::Matrix4d Matrix4x4;

// values
static PBRT_CONSTEXPR Float ShadowEpsilon = 0.0001f;
static PBRT_CONSTEXPR Float Pi = 3.14159265358979323846;
static PBRT_CONSTEXPR Float InvPi = 0.31830988618379067154;
static PBRT_CONSTEXPR Float Inv2Pi = 0.15915494309189533577;
static PBRT_CONSTEXPR Float Inv4Pi = 0.07957747154594766788;
static PBRT_CONSTEXPR Float PiOver2 = 1.57079632679489661923;
static PBRT_CONSTEXPR Float PiOver4 = 0.78539816339744830961;
static PBRT_CONSTEXPR Float Sqrt2 = 1.41421356237309504880;
static constexpr Float MaxFloat = std::numeric_limits<Float>::max();
static constexpr Float Infinity = std::numeric_limits<Float>::infinity();
#define MachineEpsilon (std::numeric_limits<Float>::epsilon() * 0.5)
