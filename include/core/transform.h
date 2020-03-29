#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "core/pbrt.h"
#include "core/bounding_boxes.h"



// this is a class of transformation matrices
// they are all 4x4 because we are in 3D world and extra dimension helps
// keep things in lienar!
class Transform
{
private:
	Eigen::Matrix4d m, mInv;


public:
	Transform() : m(Eigen::Matrix4d::Identity()), mInv(m.inverse()) {}
	Transform(const Eigen::Matrix4d &m) : m(m), mInv(m.inverse()) {}
	Transform(const Eigen::Matrix4d &m, const Eigen::Matrix4d &mInv)
		: m(m), mInv(mInv) {}

	// get the transform matrices of transpose of given transform
	friend Transform Transpose(const Transform &t);

	// transform a ray
	Ray operator()(const Ray& r) const;
		
	// TODO: p95
	Bounds3d operator()(const Bounds3d& b) const;

	// composition of transformation
	Transform operator*(const Transform& t2) const;

	// TODO: p88
	bool HasScale() const;
	
	// determine the handedness of the transformation
	// if determiant is < 0, then the handeness is changed...somehow?
	bool SwapsHandedness() const;

	SurfaceInteraction operator()(const SurfaceInteraction& si) const;



	~Transform() {}

private:
	

	

};

// give the translation matrices of given changes as delta vector
Transform Translate(const Eigen::Vector3d &delta);

// give the scale matrix of scale on x y and z
Transform Scale(double x, double y, double z);

// theta is the radian
Transform RotateX(double theta);
Transform RotateY(double theta);
Transform RotateZ(double theta);

// TODO: p90
// rotate around axis vector 
Transform Rotate(double theta, const Eigen::Vector3d &axis);

// give the transform matrix from world space to camera space
// Ex: can be used to transform camera rays to world space (with the inverse of it)
// pos, look and up in world space
// pos is the location
// look is the point that the camera is looking at
// up is the up of camera
Transform LookAt(const Eigen::Vector3d& pos, const Eigen::Vector3d& look, const Eigen::Vector3d& up);