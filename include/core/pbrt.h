// flags
#pragma once
#define _USE_MATH_DEFINES
#define PBRT_IS_WINDOWS

#ifndef PBRT_L1_CACHE_LINE_SIZE
#define PBRT_L1_CACHE_LINE_SIZE 64
#endif

// includes
// std libs
#include <math.h>
#include <iostream>
#include <vector>
#include <memory>
#include <limits>
#include <map>

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
class AreaLight;

class VisibilityTester;

class Primitive;

class Ray;
class RayDifferential;

class Transform;
class AnimatedTransform;

class Medium;
class MediumInterface;

template <typename T> class Texture;

class Shape;
class Triangle;
struct TriangleMesh;

class Interaction;
class SurfaceInteraction;

class Scene;

class Sampler;

class Camera;

class MemoryArena;

class BxDF;
class BSDF;
class BSSRDF;

class MicrofacetDistribution;
class BeckmannDistribution;
class TrowbridgeReitzDistribution;

class Material;

enum class TransportMode;

class RNG;

class Filter;

class ParamSet;
template <typename T> struct ParamSetItem;

class Integrator;


// definitions
// TODO: may have flag for float and double
typedef double Float;
// TODO: p1086
typedef Float AtomicFloat;

typedef Bounds2<Float> Bounds2f;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<Float> Bounds3f;
typedef Bounds3<int> Bounds3i;

// if change type of spectrum, need recompile
// may change later
typedef RGBSpectrum Spectrum;
//typedef SampledSpectrum Spectrum;


// using eigen for pbrt vectors/points/etc
template <typename T>
struct Vector2 : public Eigen::Matrix <T, 2, 1> 
{
	using Eigen::Matrix <T, 2, 1>::Matrix;


	template <typename U>
	explicit Vector2(const Vector2<U>& p)
	{
		this->x() = p.x();
		this->y() = p.y();
	}

};

template <typename T>
struct Point2 : public Eigen::Matrix <T, 2, 1> 
{
	using Eigen::Matrix <T, 2, 1>::Matrix;

	template <typename U>
	explicit Point2(const Point2<U>& p)
	{
		this->x() = p.x();
		this->y() = p.y();
	}

};

template <typename T>
struct Vector3 : public Eigen::Matrix <T, 3, 1> 
{
	using Eigen::Matrix <T, 3, 1>::Matrix;

	template <typename U>
	explicit Vector3(const Vector3<U>& p)
	{
		this->x() = p.x();
		this->y() = p.y();
		this->z() = p.z();
	}
};

template <typename T>
struct Point3 : public Eigen::Matrix <T, 3, 1> 
{
	using Eigen::Matrix <T, 3, 1>::Matrix;

	template <typename U>
	explicit Point3(const Point3<U>& p)
	{
		this->x() = p.x();
		this->y() = p.y();
		this->z() = p.z();
	}

};

template <typename T>
struct Normal3 : public Eigen::Matrix <T, 3, 1> 
{
	using Eigen::Matrix <T, 3, 1>::Matrix;

	template <typename U>
	explicit Normal3(const Normal3<U>& p)
	{
		this->x() = p.x();
		this->y() = p.y();
		this->z() = p.z();
	}
};

typedef Vector2<int> Vector2i;
typedef Vector2<Float> Vector2f;
typedef Point2<int> Point2i;
typedef Point2<Float> Point2f;

typedef Normal3<int> Normal3i;
typedef Normal3<Float> Normal3f;
typedef Vector3<int> Vector3i;
typedef Vector3<Float> Vector3f;
typedef Point3<int> Point3i;
typedef Point3<Float> Point3f;

typedef Eigen::Matrix4d Matrix4x4;

struct Options 
{
	int nThreads = 0;
	bool quickRender = false;
	bool quiet = false, verbose = false;
	std::string imageFile;
};

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



//#ifdef PBRT_FLOAT_IS_DOUBLE
static const Float OneMinusEpsilon = 0x1.fffffffffffffp-1;
//#else
//static const Float OneMinusEpsilon = 0x1.fffffep-1;
//#endif