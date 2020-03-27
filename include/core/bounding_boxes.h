#pragma once

#include <limits>
#include <algorithm>

#include <Eigen/Core>


template <typename T> class Bounds2
{
public:
	// vector2T
	Eigen::Matrix< T, 2, 1> pMin, pMax;

public:
	Bounds2();
	~Bounds2();

private:

};

template <typename T> class  Bounds3
{
public:
	// vector3T
	Eigen::Matrix< T, 3, 1> pMin, pMax;


public:
	Bounds3()
	{
		T minNum = std::numeric_limits<T>::lowest();
		T maxNum = std::numeric_limits<T>::max();
		pMin = Eigen::Matrix< T, 3, 1>(maxNum, maxNum, maxNum);
		pMax = Eigen::Matrix< T, 3, 1>(minNum, minNum, minNum);
	}

	Bounds3(const Eigen::Matrix< T, 3, 1>& p) : pMin(p), pMax(p) { }

	Bounds3(const Eigen::Matrix< T, 3, 1>& p1, const Eigen::Matrix< T, 3, 1>& p2)
		: pMin(std::min(p1.x(), p2.x()), std::min(p1.y(), p2.y()),
			std::min(p1.z(), p2.z())),
		pMax(std::max(p1.x(), p2.x()), std::max(p1.y(), p2.y()),
			std::max(p1.z(), p2.z())) {
	}

	~Bounds3() {}

	// select between pMin and pMax based on the value of i.
	const Eigen::Matrix< T, 3, 1>& operator[](int i) const
	{
		return i == 0 ? pMin : pMax;
	}

	Eigen::Matrix< T, 3, 1>& operator[](int i)
	{
		return i == 0 ? pMin : pMax;
	}

	Eigen::Matrix< T, 3, 1> Corner(int corner) const
	{
		return Eigen::Matrix< T, 3, 1>((*this)[(corner & 1)].x(),
			(*this)[(corner & 2) ? 1 : 0].y(),
			(*this)[(corner & 4) ? 1 : 0].z());
	}
	// not using loop like in the hw4...

	template <typename T> Bounds3<T> Union(const Bounds3<T>& b, const Eigen::Matrix< T, 3, 1>& p)
	{
		return Bounds3<T>(Eigen::Matrix< T, 3, 1>(std::min(b.pMin.x(), p.x()),
			std::min(b.pMin.y(), p.y()),
			std::min(b.pMin.z(), p.z())),
			Eigen::Matrix< T, 3, 1>(std::max(b.pMax.x(), p.x()),
				std::max(b.pMax.y(), p.y()),
				std::max(b.pMax.z(), p.z())));
	}


	template <typename T> Bounds3<T> Union(const Bounds3<T>& b1, const Bounds3<T>& b2)
	{
		return Bounds3<T>(Eigen::Matrix<T, 3, 1>(std::min(b1.pMin.x(), b2.pMin.x()),
			std::min(b1.pMin.y(), b2.pMin.y()),
			std::min(b1.pMin.z(), b2.pMin.z())),
			Eigen::Matrix< T, 3, 1>(std::max(b1.pMax.x(), b2.pMax.x()),
				std::max(b1.pMax.y(), b2.pMax.y()),
				std::max(b1.pMax.z(), b2.pMax.z())));
	}


	template <typename T> Bounds3<T> Intersect(const Bounds3<T>& b1, const Bounds3<T>& b2)
	{
		return Bounds3<T>(Eigen::Matrix<T, 3, 1>(std::max(b1.pMin.x(), b2.pMin.x()),
			std::max(b1.pMin.y() b2.pMin.y()),
			std::max(b1.pMin.z(), b2.pMin.z())),
			Eigen::Matrix<T, 3, 1>(std::min(b1.pMax.x(), b2.pMax.x()),
				std::min(b1.pMax.y(), b2.pMax.y()),
				std::min(b1.pMax.z(), b2.pMax.z())));
	}

	// TODO: maybe can simplify with Eigen
	template <typename T> bool Overlaps(const Bounds3<T>& b1, const Bounds3<T>& b2)
	{
		bool x = (b1.pMax.x() >= b2.pMin.x()) && (b1.pMin.x() <= b2.pMax.x());
		bool y = (b1.pMax.y() >= b2.pMin.y()) && (b1.pMin.y() <= b2.pMax.y());
		bool z = (b1.pMax.z() >= b2.pMin.z()) && (b1.pMin.z() <= b2.pMax.z());
		return (x && y && z);
	}

	// TODO: maybe can simplify with Eigen
	template <typename T> bool Inside(const Eigen::Matrix<T, 3, 1>& p, const Bounds3<T>& b)
	{
		return (p.x() >= b.pMin.x() && p.x() <= b.pMax.x() &&
			p.y() >= b.pMin.y() && p.y() <= b.pMax.y() &&
			p.z() >= b.pMin.z() && p.z() <= b.pMax.z());
	}


	template <typename T> bool InsideExclusive(const Eigen::Matrix<T, 3, 1>& p, const Bounds3<T>& b)
	{
		return (p.x() >= b.pMin.x() && p.x() < b.pMax.x() &&
			p.y() >= b.pMin.y() && p.y() < b.pMax.y() &&
			p.z() >= b.pMin.z() && p.z() < b.pMax.z());
	}


	template <typename T, typename U> inline Bounds3<T>	Expand(const Bounds3<T>& b, U delta)
	{
		return Bounds3<T>(b.pMin - Eigen::Matrix<T, 3, 1>(delta, delta, delta),
			b.pMax + Eigen::Matrix<T, 3, 1>(delta, delta, delta));
	}

	Eigen::Matrix<T, 3, 1> Diagonal() const { return pMax - pMin; }

	// hhat is this for??
	T SurfaceArea() const
	{
		Eigen::Matrix<T, 3, 1> d = Diagonal();
		return 2 * (d.x() * d.y() + d.x() * d.z() + d.y * d.z());
	}

	T Volume() const
	{
		Eigen::Matrix<T, 3, 1> d = Diagonal();
		return d.x() * d.y() * d.z();
	}

	int MaximumExtent() const
	{
		Eigen::Matrix<T, 3, 1> d = Diagonal();
		if (d.x() > d.y() && d.x() > d.z())
			return 0;
		else if (d.y() > d.z())
			return 1;
		else
			return 2;
	}

	// linear interoplation between conrners and given t
	Eigen::Matrix<T, 3, 1> Lerp(const Eigen::Matrix<T, 3, 1>& t) const
	{
		T x = (1 - t.x()) * pMin.x() + t.x() * pMax.x();
		T y = (1 - t.y()) * pMin.y() + t.y() * pMax.y();
		T z = (1 - t.z()) * pMin.z() + t.z() * pMax.z();
		return Eigen::Matrix<T, 3, 1>(x, y, z);
	}

	//TODO: ...

private:

};

typedef Bounds2<double> Bounds2d;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<double> Bounds3d;
typedef Bounds3<int> Bounds3i;


