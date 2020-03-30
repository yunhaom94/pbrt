#include "core/triangle_mesh.h"
#include "core/transform.h"

TriangleMesh::TriangleMesh(const Transform& ObjectToWorld, 
	int nTriangles, 
	const int& vertexIndices, 
	int nVertices, 
	const Eigen::Vector3d* P, 
	const Eigen::Vector3d* S, 
	const Eigen::Vector3d* N, 
	const Eigen::Vector2d* UV,
	const std::shared_ptr<Texture<double>>& alphaMask) :
	nTriangles(nTriangles), 
	nVertices(nVertices),
	vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles),
	alphaMask(alphaMask)
{

	// Transform mesh vertices to world space
	p.reset(new Eigen::Vector3d[nVertices]);
	for (int i = 0; i < nVertices; ++i)
		p[i] = ObjectToWorld(P[i]);
}
