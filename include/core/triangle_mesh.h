
#pragma once

#include "core/pbrt.h"


// V ~= 2F for closed meshes
struct TriangleMesh
{
	const int nTriangles, nVertices;
	// F, where n * i, n * i + 1, n * i + 2 is a triangle
	std::vector<int> vertexIndices;
	// V as an array of 3d vector
	std::unique_ptr<Eigen::Vector3d[]> p;
	// N as an array of 3d vector
	std::unique_ptr<Eigen::Vector3d[]> n;
	// tagent vector as an array of 3d vector
	std::unique_ptr<Eigen::Vector3d[]> s;
	// UV as an array of 2d vector
	std::unique_ptr<Eigen::Vector2d[]> uv;
	std::shared_ptr<Texture<double>> alphaMask;


	TriangleMesh(const Transform& ObjectToWorld,
		int nTriangles, const int& vertexIndices, int nVertices,
		const Eigen::Vector3d* P, const  Eigen::Vector3d* S, const  Eigen::Vector3d* N,
		const  Eigen::Vector2d* UV,
		const std::shared_ptr<Texture<double>>& alphaMask);
};