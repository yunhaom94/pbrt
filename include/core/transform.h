#pragma once

#include "core/pbrt.h"
#include "utlis/utlis.h"

// this is a class of transformation matrices
// they are all 4x4 because we are in 3D world and extra dimension helps
// keep things in lienar!
class Transform
{
private:
	Matrix4x4 m, mInv;


public:
	Transform() : m(Matrix4x4::Identity()), mInv(m.inverse()) {}
	Transform(const Float mat[4][4]) 
	{
		m << mat[0][0], mat[0][1], mat[0][2], mat[0][3],
			mat[1][0], mat[1][1], mat[1][2], mat[1][3],
			mat[2][0], mat[2][1], mat[2][2], mat[2][3],
			mat[3][0], mat[3][1], mat[3][2], mat[3][3];
		mInv = (m.inverse());
	}
	Transform(const Matrix4x4 &m) : m(m), mInv(m.inverse()) {}
	Transform(const Matrix4x4 &m, const Matrix4x4 &mInv)
		: m(m), mInv(mInv) {}

	~Transform() {}

	// return a new transform with matrices inversed
	friend Transform Inverse(const Transform& t);

	// return a new transform with matrices of transpose of given transform
	friend Transform Transpose(const Transform &t);

	// transform a point 
	template <typename T>
	inline Point3<T> operator()(const Point3<T>& p) const;

	template <typename T>
	inline Point3<T> operator()(const Point3<T>& p, Vector3<T>* pError) const;

	template<typename T>
	inline Point3<T> operator()(const Point3<T>& pt, const Vector3<T>& ptError, Vector3<T>* absError) const;

	template <typename T>
	inline Vector3<T> operator()(const Vector3<T>& v) const;

	template <typename T>
	inline Vector3<T> operator()(const Vector3<T>& v, Vector3<T>* vTransError) const;

	template <typename T>
	inline Vector3<T> operator()(const Vector3<T>& v, const Vector3<T>& vError, Vector3<T>* vTransError) const;

	// transform a normal
	template <typename T>
	inline Normal3<T> operator()(const Normal3<T>& n) const;
		
	// transform a ray
	Ray operator()(const Ray& r) const;

	Ray operator()(const Ray& r, Vector3f* oError,
		Vector3f* dError) const;

	Ray operator()(const Ray& r, const Vector3f& oErrorIn,
		const Vector3f& dErrorIn, Vector3f* oErrorOut,
		Vector3f* dErrorOut) const;

	RayDifferential operator()(const RayDifferential& r) const;

		
	// transform a bound box
	Bounds3f operator()(const Bounds3f& b) const;

	// transform a surface interaction
	SurfaceInteraction operator()(const SurfaceInteraction& si) const;

	// composition of transformation
	Transform operator*(const Transform& t2) const;

	// return if a transformation has a scaling term in it
	// by transform 3 basis vectors and see if any of their 
	// lengths are appreciably different from one.
	bool HasScale() const;
	
	// determine the handedness of the transformation
	// if determiant is < 0, then the handeness is changed...somehow?
	bool SwapsHandedness() const;

	bool IsIdentity() const
	{
		return m.isIdentity();
	}
};

// TODO: p97
class AnimatedTransform
{
public:
	AnimatedTransform() {}
	~AnimatedTransform() {}

	void Interpolate(Float time, Transform* t) const {}

	Ray operator()(const Ray& r) const;
	RayDifferential operator()(const RayDifferential& r) const;
	Point3f operator()(Float time, const Point3f& p) const;
	Vector3f operator()(Float time, const Vector3f& v) const;

private:

};


// give the translation matrices of given changes as delta vector
Transform Translate(const Vector3f &delta);

// give the scale matrix of scale on x y and z
Transform Scale(Float x, Float y, Float z);

// theta is the radian
Transform RotateX(Float theta);
Transform RotateY(Float theta);
Transform RotateZ(Float theta);

// rotate around axis vector 
// theta is in radians
Transform Rotate(Float theta, const Vector3f& axis);

// give the transform matrix from world space to camera space
// Ex: can be used to transform camera rays to world space (with the inverse of it)
// pos, look and up in world space
// pos is the location
// look is the point that the camera is looking at
// up is the up of camera
Transform LookAt(const Vector3f& pos, const Vector3f& look, const Vector3f& up);

Transform Orthographic(Float zNear, Float zFar);

Transform Perspective(Float fov, Float n, Float f);

// TODO: p1124
inline void TransformCache() {}

template<typename T>
inline Point3<T> Transform::operator()(const Point3<T>& p) const
{
	
	T x = p.x(), y = p.y(), z = p.z();
	T wp = m.row(3)[0] * x + m.row(3)[1] * y + m.row(3)[2] * z + m.row(3)[3];

	Point3<T> ret = (Eigen::Matrix <T, 4, 1>(x, y, z, 0).transpose() * m).head(3);

	if (wp == 1)
		return ret;
	else 
		return ret / wp;
	
}

template <typename T>
inline Point3<T> Transform::operator()(const Point3<T>& p,
	Vector3<T>* pError) const
{
	T x = p.x(), y = p.y(), z = p.z();
	// Compute transformed coordinates from point _pt_
	T xp = (m.row(0)[0] * x + m.row(0)[1] * y) + (m.row(0)[2] * z + m.row(0)[3]);
	T yp = (m.row(1)[0] * x + m.row(1)[1] * y) + (m.row(1)[2] * z + m.row(1)[3]);
	T zp = (m.row(2)[0] * x + m.row(2)[1] * y) + (m.row(2)[2] * z + m.row(2)[3]);
	T wp = (m.row(3)[0] * x + m.row(3)[1] * y) + (m.row(3)[2] * z + m.row(3)[3]);

	// Compute absolute error for transformed point
	T xAbsSum = (std::abs(m.row(0)[0] * x) + std::abs(m.row(0)[1] * y) +
				 std::abs(m.row(0)[2] * z) + std::abs(m.row(0)[3]));
	T yAbsSum = (std::abs(m.row(1)[0] * x) + std::abs(m.row(1)[1] * y) +
				 std::abs(m.row(1)[2] * z) + std::abs(m.row(1)[3]));
	T zAbsSum = (std::abs(m.row(2)[0] * x) + std::abs(m.row(2)[1] * y) +
		         std::abs(m.row(2)[2] * z) + std::abs(m.row(2)[3]));
	*pError = gamma(3) * Eigen::Matrix<T, 3, 1>(xAbsSum, yAbsSum, zAbsSum);
	//CHECK_NE(wp, 0);
	if (wp == 1)
		return Point3<T>(xp, yp, zp);
	else
		return Point3<T>(xp, yp, zp) / wp;
}

template <typename T>
inline Point3<T> Transform::operator()(const Point3<T>& pt,
	const Vector3<T>& ptError,
	Vector3<T>* absError) const
{
	T x = pt.x(), y = pt.y(), z = pt.z();
	T xp = (m.row(0)[0] * x + m.row(0)[1] * y) + (m.row(0)[2] * z + m.row(0)[3]);
	T yp = (m.row(1)[0] * x + m.row(1)[1] * y) + (m.row(1)[2] * z + m.row(1)[3]);
	T zp = (m.row(2)[0] * x + m.row(2)[1] * y) + (m.row(2)[2] * z + m.row(2)[3]);
	T wp = (m.row(3)[0] * x + m.row(3)[1] * y) + (m.row(3)[2] * z + m.row(3)[3]);
	absError->x() =
		(gamma(3) + (T)1) *
		(std::abs(m.row(0)[0]) * ptError.x() + std::abs(m.row(0)[1]) * ptError.y() +
			std::abs(m.row(0)[2]) * ptError.z()) +
		gamma(3) * (std::abs(m.row(0)[0] * x) + std::abs(m.row(0)[1] * y) +
			std::abs(m.row(0)[2] * z) + std::abs(m.row(0)[3]));
	absError->y() =
		(gamma(3) + (T)1) *
		(std::abs(m.row(1)[0]) * ptError.x() + std::abs(m.row(1)[1]) * ptError.y() +
			std::abs(m.row(1)[2]) * ptError.z()) +
		gamma(3) * (std::abs(m.row(1)[0] * x) + std::abs(m.row(1)[1] * y) +
			std::abs(m.row(1)[2] * z) + std::abs(m.row(1)[3]));
	absError->z() =
		(gamma(3) + (T)1) *
		(std::abs(m.row(2)[0]) * ptError.x() + std::abs(m.row(2)[1]) * ptError.y() +
			std::abs(m.row(2)[2]) * ptError.z()) +
		gamma(3) * (std::abs(m.row(2)[0] * x) + std::abs(m.row(2)[1] * y) +
			std::abs(m.row(2)[2] * z) + std::abs(m.row(2)[3]));
	//CHECK_NE(wp, 0);
	if (wp == 1)
		return Point3<T>(xp, yp, zp);
	else
		return Point3<T>(xp, yp, zp) / wp;
}

template<typename T>
inline Vector3<T> Transform::operator()(const Vector3<T>& v) const
{
	T x = v.x(), y = v.y(), z = v.z();
	Vector3<T> ret = (Eigen::Matrix <T, 4, 1>(x, y, z, 0).transpose() * m).head(3);
	return ret;
}

template<typename T>
inline Vector3<T> Transform::operator()(const Vector3<T>& v, Vector3<T>* absError) const
{
	T x = v.x(), y = v.y(), z = v.z();
	absError->x() =
		gamma(3) * (std::abs(m.row(0)[0] * v.x()) + std::abs(m.row(0)[1] * v.y()) +
			std::abs(m.row(0)[2] * v.z()));
	absError->y() =
		gamma(3) * (std::abs(m.row(1)[0] * v.x()) + std::abs(m.row(1)[1] * v.y()) +
			std::abs(m.row(1)[2] * v.z()));
	absError->z() =
		gamma(3) * (std::abs(m.row(2)[0] * v.x()) + std::abs(m.row(2)[1] * v.y()) +
			std::abs(m.row(2)[2] * v.z()));
	return Vector3<T>(m.row(0)[0] * x + m.row(0)[1] * y + m.row(0)[2] * z,
		m.row(1)[0] * x + m.row(1)[1] * y + m.row(1)[2] * z,
		m.row(2)[0] * x + m.row(2)[1] * y + m.row(2)[2] * z);
}

template<typename T>
inline Vector3<T> Transform::operator()(const Vector3<T>& v, const Vector3<T>& vError, Vector3<T>* absError) const
{
	T x = v.x(), y = v.y(), z = v.z();
	absError->x() =
		(gamma(3) + (T)1) *
		(std::abs(m.row(0)[0]) * vError.x() + std::abs(m.row(0)[1]) * vError.y() +
			std::abs(m.row(0)[2]) * vError.z()) +
		gamma(3) * (std::abs(m.row(0)[0] * v.x()) + std::abs(m.row(0)[1] * v.y()) +
			std::abs(m.row(0)[2] * v.z()));
	absError->y() =
		(gamma(3) + (T)1) *
		(std::abs(m.row(1)[0]) * vError.x() + std::abs(m.row(1)[1]) * vError.y() +
			std::abs(m.row(1)[2]) * vError.z()) +
		gamma(3) * (std::abs(m.row(1)[0] * v.x()) + std::abs(m.row(1)[1] * v.y()) +
			std::abs(m.row(1)[2] * v.z()));
	absError->z() =
		(gamma(3) + (T)1) *
		(std::abs(m.row(2)[0]) * vError.x() + std::abs(m.row(2)[1]) * vError.y() +
			std::abs(m.row(2)[2]) * vError.z()) +
		gamma(3) * (std::abs(m.row(2)[0] * v.x()) + std::abs(m.row(2)[1] * v.y()) +
			std::abs(m.row(2)[2] * v.z()));
	return Vector3<T>(m.row(0)[0] * x + m.row(0)[1] * y + m.row(0)[2] * z,
		m.row(1)[0] * x + m.row(1)[1] * y + m.row(1)[2] * z,
		m.row(2)[0] * x + m.row(2)[1] * y + m.row(2)[2] * z);
}

template<typename T>
inline Normal3<T> Transform::operator()(const Normal3<T>&n) const
{
	T x = n.x(), y = n.y(), z = n.z();
	Normal3<T> ret = (Eigen::Matrix <T, 4, 1>(x, y, z, 0).transpose() * mInv).head(3);
	return ret;
}



