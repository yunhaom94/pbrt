#pragma once

#include "core/pbrt.h"

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
	inline Eigen::Matrix <T, 3, 1> TransformPoint(const Eigen::Matrix <T, 3, 1>& p) const;

	// transform a vector
	template <typename T>
	inline Eigen::Matrix <T, 3, 1> TransformVector(const Eigen::Matrix <T, 3, 1>& p) const;

	// transform a normal
	template <typename T>
	inline Eigen::Matrix <T, 3, 1> TransformNormal(const Eigen::Matrix <T, 3, 1>& p) const;
		
	// transform a ray
	Ray operator()(const Ray& r) const;
		
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
};

// TODO: p97
class AnimatedTransform
{
public:
	AnimatedTransform() {}
	~AnimatedTransform() {}

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

// TODO: p1124
inline void TransformCache() {}

template<typename T>
inline Eigen::Matrix<T, 3, 1> Transform::TransformPoint(const Eigen::Matrix<T, 3, 1>& p) const
{
	
	T x = p.x(), y = p.y(), z = p.z();
	T wp = m.row(3)[0] * x + m.row(3)[1] * y + m.row(3)[2] * z + m.row(3)[3];

	Eigen::Matrix <T, 3, 1> ret = (Eigen::Matrix <T, 4, 1>(x, y, z, 0).transpose() * m).head(3);

	if (wp == 1)
		return ret;
	else 
		return ret / wp;
	
}

template<typename T>
inline Eigen::Matrix<T, 3, 1> Transform::TransformVector(const Eigen::Matrix<T, 3, 1>& p) const
{
	T x = p.x(), y = p.y(), z = p.z();
	Eigen::Matrix <T, 3, 1> ret = (Eigen::Matrix <T, 4, 1>(x, y, z, 0).transpose() * m).head(3);
}

template<typename T>
inline Eigen::Matrix<T, 3, 1> Transform::TransformNormal(const Eigen::Matrix<T, 3, 1>& p) const
{
	T x = p.x(), y = p.y(), z = p.z();
	Eigen::Matrix <T, 3, 1> ret = (Eigen::Matrix <T, 4, 1>(x, y, z, 0).transpose() * mInv).head(3);
}

