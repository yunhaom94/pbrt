#pragma once

#include "core/pbrt.h"
#include "core/shape.h"


class Triangle : public Shape
{
public:
	const std::shared_ptr<TriangleMesh> mesh;
	
	// pointer into the vertex list
	// TODO: may change later to actual vertices
	const int* v;

public:
	Triangle(const Transform* ObjectToWorld, const Transform* WorldToObject,
		bool reverseOrientation,
		const std::shared_ptr<TriangleMesh>& mesh, int triNumber);

	//TODO: p156
	std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
		const Transform* ObjectToWorld, 
		const Transform* WorldToObject,
		bool reverseOrientation, 
		int nTriangles,
		const int* vertexIndices,
		int nVertices, 
		const Eigen::Vector3d* p,
		const Eigen::Vector3d* s, 
		const Eigen::Vector3d* n, 
		const Eigen::Vector2d* uv,
		const std::shared_ptr<Texture<double>>& alphaMask);

	// create a bounding box around the triangle
	Bounds3d ObjectBound() const;

	// calcuate intersect between triangle and a ray
	// used algorithm in hw4
	bool Triangle::Intersect(const Ray& ray, double& tHit,
		SurfaceInteraction* isect, bool testAlphaTexture) const;

	~Triangle() {};

private:

};

