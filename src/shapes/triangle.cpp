#include "shapes/triangle.h"

#include "core/triangle_mesh.h"
#include "core/transform.h"
#include "core/bounding_boxes.h"
#include "core/interaction.h"

/*
Triangle::Triangle(const Transform* ObjectToWorld, 
	const Transform* WorldToObject, 
	bool reverseOrientation, 
	const std::shared_ptr<TriangleMesh>& mesh, int triNumber) : 
	Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh)
{
	v = &mesh->vertexIndices[3 * triNumber];
}

std::vector<std::shared_ptr<Shape>> Triangle::CreateTriangleMesh(const Transform* ObjectToWorld, const Transform* WorldToObject,
	bool reverseOrientation, int nTriangles, 
	const int* vertexIndices, int nVertices, 
	const Eigen::Vector3d* p, const Eigen::Vector3d* s, 
	const Eigen::Vector3d* n, const Eigen::Vector2d* uv, 
	const std::shared_ptr<Texture<double>>& alphaMask)
{
	return std::vector<std::shared_ptr<Shape>>();
}

Bounds3f Triangle::ObjectBound() const
{
	const Eigen::Vector3d& p0 = mesh->p[v[0]];
	const Eigen::Vector3d& p1 = mesh->p[v[1]];
	const Eigen::Vector3d& p2 = mesh->p[v[2]];
    return Bounds3f();
	//return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),
		//(*WorldToObject)(p2));
}

bool Triangle::Intersect(const Ray& ray, double& tHit, SurfaceInteraction* isect, bool testAlphaTexture) const
{
    double min_t = 0;

    Point3f a = (*ObjectToWorld)(mesh->p[v[0]]);
    Point3f b = (*ObjectToWorld)(mesh->p[v[1]]);
    Point3f c = (*ObjectToWorld)(mesh->p[v[2]]);

    Eigen::Vector3d e = ray.o;
    Eigen::Vector3d d = ray.d;

    Eigen::Vector3d t1 = a - b;
    Eigen::Vector3d t2 = a - c;

    Eigen::Vector3d f = a - e;

    // book 4.4
    double m = t1[0] * (t2[1] * d[2] - d[1] * t2[2]) +
        t1[1] * (d[0] * t2[2] - t2[0] * d[2]) +
        t1[2] * (t2[0] * d[1] - t2[1] * d[0]);

    double t_res = -((t2[2] * (t1[0] * f[1] - f[0] * t1[1]) +
        t2[1] * (f[0] * t1[2] - t1[0] * f[2]) +
        t2[0] * (t1[1] * f[2] - f[1] * t1[2])) /
        m);

    if (t_res < min_t || t_res > ray.tMax)
        return false;

    double gamma = (d[2] * (t1[0] * f[1] - f[0] * t1[1]) +
        d[1] * (f[0] * t1[2] - t1[0] * f[2]) +
        d[0] * (t1[1] * f[2] - f[1] * t1[2])) / m;

    if (gamma < 0 || gamma > 1)
        return false;

    double beta = (f[0] * (t2[1] * d[2] - d[1] * t2[2]) +
        f[1] * (d[0] * t2[2] - t2[0] * d[2]) +
        f[2] * (t2[0] * d[1] - t2[1] * d[0])) /
        m;

    if (beta < 0 || beta > 1 - gamma)
        return false;


    Eigen::Vector3d ba = b - a;
    Eigen::Vector3d ca = c - a;

    tHit = t_res;
    Eigen::Vector3d n = ba.cross(ca);

    // derivatives
    Eigen::Vector3d dpdu(0, 0, 0), dpdv(0, 0, 0), pError(0,0,0);
    Eigen::Vector2d uvHit(0,0);
   // TODO: p157 NOT FINISHED

    Eigen::Vector3d pHit = e + t_res * d;

    *isect = SurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv,
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), ray.time, this);
    isect->n = isect->shading.n = n;



	return true;
}
*/