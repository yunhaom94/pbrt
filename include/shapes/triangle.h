#pragma once

#include "core/pbrt.h"
#include "core/shape.h"

// TODO: remove this
#include "core/interaction.h"
template <typename T> class Texture
{
public:
	Texture() {}
	~Texture() {}

	int Evaluate(SurfaceInteraction s) { return 0; }

};

// V ~= 2F for closed meshes
struct TriangleMesh
{
	// # of triangles and vertices
	const int nTriangles, nVertices;
	// F, where n * i, n * i + 1, n * i + 2 is a triangle
	std::vector<int> vertexIndices;
	// V as an array of 3d vector
	std::unique_ptr<Point3f[]> p;
	// N as an array of 3d vector
	std::unique_ptr<Normal3f[]> n;
	// tagent vector as an array of 3d vector
	std::unique_ptr<Vector3f[]> s;
	// UV as an array of 2d vector
	std::unique_ptr<Point2f[]> uv;
	std::shared_ptr<Texture<Float>> alphaMask;


	TriangleMesh(const Transform& ObjectToWorld,
		int nTriangles, const int* vertexIndices, int nVertices,
		const Point3f* P, const Vector3f* S, const Normal3f* N,
		const Point2f* UV,
		const std::shared_ptr<Texture<Float>>& alphaMask);

	~TriangleMesh() {}
};

class Triangle : public Shape
{
public:
	const std::shared_ptr<TriangleMesh> mesh;
	
	// a pointer to the first vertex index
	// get all three with v[0], v[1], v[2] 
	const int* v;

public:
	Triangle(const Transform* ObjectToWorld, const Transform* WorldToObject,
		bool reverseOrientation,
		const std::shared_ptr<TriangleMesh>& mesh, int triNumber);

	// create a bunch of triangles and init a mesh and add 
	// reference to the triangles
	std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
		const Transform* ObjectToWorld, 
		const Transform* WorldToObject,
		bool reverseOrientation, 
		int nTriangles,
		const int* vertexIndices,
		int nVertices,
		const Point3f* p,
		const Vector3f* s,
		const Normal3f* n, 
		const Point2f* uv,
		const std::shared_ptr<Texture<Float>>& alphaMask);

	~Triangle() {};

	// create a bounding box around the triangle
	Bounds3f ObjectBound() const;

	// calculate intersect between triangle and a ray
	// we use ray coordinate space this time
	bool Intersect(const Ray& ray, Float *tHit,
		SurfaceInteraction* isect, bool testAlphaTexture) const;

	// TODO: maybe later used algorithm in hw4
	//bool Intersect418(const Ray& ray, double& tHit,
	//	SurfaceInteraction* isect, bool testAlphaTexture) const;

	bool IntersectP(const Ray& r, bool testAlphaTexture) const;

	void GetUVs(Point2f uv[3]) const;

	double Area() const;

	Interaction Triangle::Sample(const Point2f& u) const;

private:

};

