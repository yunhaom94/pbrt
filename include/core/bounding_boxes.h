#pragma once

#include "core/pbrt.h"
#include "core/ray.h"

// DO NOT FORMAT THIS DOCUMENT IN VS

// TODO:
template <typename T> class Bounds2
{
public:
	// vector2T
	Eigen::Matrix< T, 2, 1> pMin, pMax;

public:
	Bounds2() {}
	Bounds2(const Eigen::Matrix< T, 2, 1>& p1, const Eigen::Matrix< T, 2, 1>& p2) : pMin(p1), pMax(p2) { }
	~Bounds2() {}


	Vector2i Diagonal()
	{
		return Vector2i(0, 0);
	}

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
		// default constructors create an empty box by 
		// setting the extent to an invalid configuration
		T minNum = std::numeric_limits<T>::lowest();
		T maxNum = std::numeric_limits<T>::max();
		pMin = Eigen::Matrix< T, 3, 1>(maxNum, maxNum, maxNum);
		pMax = Eigen::Matrix< T, 3, 1>(minNum, minNum, minNum);
	}

	Bounds3(const Eigen::Matrix< T, 3, 1>& p) : pMin(p), pMax(p) { }

	Bounds3(const Eigen::Matrix< T, 3, 1>& p1, 
		const Eigen::Matrix< T, 3, 1>& p2)
		: pMin(std::min(p1.x(), p2.x()), 
			   std::min(p1.y(), p2.y()),
			   std::min(p1.z(), p2.z())),
		  pMax(std::max(p1.x(), p2.x()), 
			   std::max(p1.y(), p2.y()),
			   std::max(p1.z(), p2.z())) { }

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

	// returns the coordinates of one of the eight corners of the BB
	Eigen::Matrix< T, 3, 1> Corner(int corner) const
	{
		return Eigen::Matrix< T, 3, 1>((*this)[(corner & 1)].x(),
			(*this)[(corner & 2) ? 1 : 0].y(),
			(*this)[(corner & 4) ? 1 : 0].z());
	}

	Eigen::Matrix<T, 3, 1> Diagonal() const { return pMax - pMin; }

	// k?
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
	// returns the index of which of the three axes is longest
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

	// linear interpolation between corners and given t
	Eigen::Matrix<T, 3, 1> Lerp(const Eigen::Matrix<T, 3, 1>& t) const
	{
		T x = (1 - t.x()) * pMin.x() + t.x() * pMax.x();
		T y = (1 - t.y()) * pMin.y() + t.y() * pMax.y();
		T z = (1 - t.z()) * pMin.z() + t.z() * pMax.z();
		return Eigen::Matrix<T, 3, 1>(x, y, z);
	}

	// Offset() returns the continuous position of a point relative to the corners of the box
	// starting at the min corner as origin.
	Eigen::Matrix<T, 3, 1> Offset(const Eigen::Matrix<T, 3, 1>& p) const 
	{
		Eigen::Matrix<T, 3, 1> o = p - pMin;
		if (pMax.x() > pMin.x()) 
			o.x() /= pMax.x() - pMin.x();
		if (pMax.y() > pMin.y()) 
			o.y() /= pMax.y() - pMin.y();
		if (pMax.z() > pMin.z()) 
			o.z() /= pMax.z() - pMin.z();
		return o;
	}

	// write a loosely bounding sphere to provided
	// center and radius
	void BoundingSphere(Point3f* center, Float* radius) const 
	{
		*center = (pMin + pMax) / 2;
		*radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
	}

	// TODO: iterator for something? p81

	// From hw4 ray_intersect_box
	template <typename T>
	bool IntersectP(const Ray& ray, double& hitt0,
		double& hitt1) const
	{
		Eigen::Vector3d o = ray.o;
		Eigen::Vector3d d = ray.d;

		T xmin = pMin[0];
		T ymin = pMin[1];
		T zmin = pMin[2];
		T xmax = pMax[0];
		T ymax = pMax[1];
		T zmax = pMax[2];

		// inside box
		if (o[0] >= xmin &&
			o[1] >= ymin &&
			o[2] >= zmin &&
			o[0] <= xmax &&
			o[1] <= ymax &&
			o[2] <= zmax)
		{
			hitt0 = 0;
			return true;
		}

		T txmin, txmax, tymin, tymax, tzmin, tzmax;

		double xrec = 1 / d[0];
		if (xrec >= 0)
		{
			txmin = xrec * (xmin - o[0]);
			txmax = xrec * (xmax - o[0]);
		}
		else
		{
			txmin = xrec * (xmax - o[0]);
			txmax = xrec * (xmin - o[0]);
		}

		double yrec = 1 / d[1];
		if (yrec >= 0)
		{
			tymin = yrec * (ymin - o[1]);
			tymax = yrec * (ymax - o[1]);
		}
		else
		{
			tymin = yrec * (ymax - o[1]);
			tymax = yrec * (ymin - o[1]);
		}

		double zrec = 1 / d[2];
		if (zrec >= 0)
		{
			tzmin = zrec * (zmin - o[2]);
			tzmax = zrec * (zmax - o[2]);
		}
		else
		{
			tzmin = zrec * (zmax - o[2]);
			tzmax = zrec * (zmin - o[2]);
		}

		double tmin = std::max(std::max(txmin, tymin), tzmin);
		double tmax = std::min(std::min(txmax, tymax), tzmax);

		if (tmin > tmax || tmin > ray.tMax)
			return false;
		else
		{
			hitt0 = tmin;
			hitt1 = tmax;
			return true;
		}
	}

	// same as last one but with recs already calculated
	template <typename T>
	bool IntersectP(const Ray& ray, const Eigen::Vector3d& invDir,
		const int dirIsNeg[3]) const
	{
		//TODO: p128
		return false;
	}

	//TODO: 


private:

};


// utility functions of bounding boxes

template <typename T> Bounds3<T> 
Union(const Bounds3<T>& b, const Eigen::Matrix< T, 3, 1>& p)
{
	return Bounds3<T>(Eigen::Matrix< T, 3, 1>(std::min(b.pMin.x(), p.x()),
											  std::min(b.pMin.y(), p.y()),
											  std::min(b.pMin.z(), p.z())),
					  Eigen::Matrix< T, 3, 1>(std::max(b.pMax.x(), p.x()),
											  std::max(b.pMax.y(), p.y()),
											  std::max(b.pMax.z(), p.z())));
}


template <typename T> Bounds3<T> 
Union(const Bounds3<T>& b1, const Bounds3<T>& b2)
{
	return Bounds3<T>(Eigen::Matrix<T, 3, 1>(std::min(b1.pMin.x(), b2.pMin.x()),
											 std::min(b1.pMin.y(), b2.pMin.y()),
											 std::min(b1.pMin.z(), b2.pMin.z())),
					  Eigen::Matrix<T, 3, 1>(std::max(b1.pMax.x(), b2.pMax.x()),
											 std::max(b1.pMax.y(), b2.pMax.y()),
											 std::max(b1.pMax.z(), b2.pMax.z())));
}


template <typename T> Bounds3<T>
Intersect(const Bounds3<T>& b1, const Bounds3<T>& b2)
{
	return Bounds3<T>(Eigen::Matrix<T, 3, 1>(std::max(b1.pMin.x(), b2.pMin.x()),
											 std::max(b1.pMin.y(), b2.pMin.y()),
											 std::max(b1.pMin.z(), b2.pMin.z())),
					  Eigen::Matrix<T, 3, 1>(std::min(b1.pMax.x(), b2.pMax.x()),
											 std::min(b1.pMax.y(), b2.pMax.y()),
											 std::min(b1.pMax.z(), b2.pMax.z())));
}

// TODO: maybe can simplify with Eigen
template <typename T> 
bool Overlaps(const Bounds3<T>& b1, const Bounds3<T>& b2)
{
	bool x = (b1.pMax.x() >= b2.pMin.x()) && (b1.pMin.x() <= b2.pMax.x());
	bool y = (b1.pMax.y() >= b2.pMin.y()) && (b1.pMin.y() <= b2.pMax.y());
	bool z = (b1.pMax.z() >= b2.pMin.z()) && (b1.pMin.z() <= b2.pMax.z());
	return (x && y && z);
}

// TODO: maybe can simplify with Eigen
template <typename T> 
bool Inside(const Eigen::Matrix<T, 3, 1>& p, const Bounds3<T>& b)
{
	return (p.x() >= b.pMin.x() && p.x() <= b.pMax.x() &&
		    p.y() >= b.pMin.y() && p.y() <= b.pMax.y() &&
		    p.z() >= b.pMin.z() && p.z() <= b.pMax.z());
}


template <typename T> 
bool InsideExclusive(const Eigen::Matrix<T, 3, 1>& p, const Bounds3<T>& b)
{
	return (p.x() >= b.pMin.x() && p.x() < b.pMax.x() &&
		    p.y() >= b.pMin.y() && p.y() < b.pMax.y() &&
		    p.z() >= b.pMin.z() && p.z() < b.pMax.z());
}


template <typename T, typename U> 
inline Bounds3<T> Expand(const Bounds3<T>& b, U delta)
{
	return Bounds3<T>(b.pMin - Eigen::Matrix<T, 3, 1>(delta, delta, delta),
		              b.pMax + Eigen::Matrix<T, 3, 1>(delta, delta, delta));
}