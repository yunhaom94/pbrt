#include "shapes/triangle.h"
#include "core/transform.h"
#include "core/bounding_boxes.h"
#include "core/interaction.h"
#include "utlis/geometry.h"



Triangle::Triangle(const Transform* ObjectToWorld, 
	const Transform* WorldToObject, 
	bool reverseOrientation, 
	const std::shared_ptr<TriangleMesh>& mesh, int triNumber) : 
	Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh)
{
	v = &mesh->vertexIndices[3 * triNumber];
}


std::vector<std::shared_ptr<Shape>> Triangle::CreateTriangleMesh(
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
    const std::shared_ptr<Texture<Float>>& alphaMask)
{
    std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
        *ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv,
        alphaMask);

    // list of triangles
    std::vector<std::shared_ptr<Shape>> tris;

    for (int i = 0; i < nTriangles; ++i)
        tris.push_back(std::make_shared<Triangle>(ObjectToWorld,
            WorldToObject, reverseOrientation, mesh, i));

    return tris;
}



Bounds3f Triangle::ObjectBound() const
{
	const Point3f& p0 = mesh->p[v[0]];
	const Point3f& p1 = mesh->p[v[1]];
	const Point3f& p2 = mesh->p[v[2]];
    return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),
        (*WorldToObject)(p2));
}

bool Triangle::Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect, bool testAlphaTexture) const
{
    // 3 vertices
    const Point3f& p0 = mesh->p[v[0]];
    const Point3f& p1 = mesh->p[v[1]];
    const Point3f& p2 = mesh->p[v[2]];

    // create a to-ray-space transformation
    // using weird permutation thing
    // with ray starts at 0,0,0
    Point3f p0t = p0 - Vector3f(ray.o);
    Point3f p1t = p1 - Vector3f(ray.o);
    Point3f p2t = p2 - Vector3f(ray.o);

    int kz = MaxDimension<Float>(ray.d.cwiseAbs());
    int kx = kz + 1; if (kx == 3) kx = 0;
    int ky = kx + 1; if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.d, kx, ky, kz);
    p0t = Permute<Float>(p0t, kx, ky, kz);
    p1t = Permute<Float>(p1t, kx, ky, kz);
    p2t = Permute<Float>(p2t, kx, ky, kz);

    Float Sx = -d.x() / d.z();
    Float Sy = -d.y() / d.z();
    Float Sz = 1.0 / d.z();
    p0t.x() += Sx * p0t.z();
    p0t.y() += Sy * p0t.z();
    p1t.x() += Sx * p1t.z();
    p1t.y() += Sy * p1t.z();
    p2t.x() += Sx * p2t.z();
    p2t.y() += Sy * p2t.z();

    Float e0 = p1t.x() * p2t.y() - p1t.y() * p2t.x();
    Float e1 = p2t.x() * p0t.y() - p2t.y() * p0t.x();
    Float e2 = p0t.x() * p1t.y() - p0t.y() * p1t.x();

    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    Float det = e0 + e1 + e2;
    if (det == 0)
        return false;

    p0t.z() *= Sz;
    p1t.z() *= Sz;
    p2t.z() *= Sz;
    Float tScaled = e0 * p0t.z() + e1 * p1t.z() + e2 * p2t.z();
    if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
        return false;

    Float invDet = 1 / det;
    Float b0 = e0 * invDet;
    Float b1 = e1 * invDet;
    Float b2 = e2 * invDet;
    Float t = tScaled * invDet;

    // calculate interaction params
    Vector3f dpdu, dpdv;
    Point2f uv[3];
    GetUVs(uv);

    Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
    Vector3f dp02 = p0 - p2, dp12 = p1 - p2;

    Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
    if (determinant == 0) 
    {
        CoordinateSystem<Float>((p2 - p0).cross(p1 - p0).normalized(), &dpdu, &dpdv);
    }
    else 
    {
        Float invdet = 1 / determinant;
        dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
        dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }

    Float xAbsSum =
        (std::abs(b0 * p0.x()) + std::abs(b1 * p1.x()) + std::abs(b2 * p2.x()));
    Float yAbsSum =
        (std::abs(b0 * p0.y()) + std::abs(b1 * p1.y()) + std::abs(b2 * p2.y()));
    Float zAbsSum =
        (std::abs(b0 * p0.z()) + std::abs(b1 * p1.z()) + std::abs(b2 * p2.z()));
    Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);


    Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
    Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

    // if hits transparent texture, just pretend no hit
    if (testAlphaTexture && mesh->alphaMask)
    {
        SurfaceInteraction isectLocal(pHit, Vector3f(0, 0, 0), uvHit,
            Vector3f(0, 0, 0), dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0),
            ray.time, this);
        if (mesh->alphaMask->Evaluate(isectLocal) == 0)
            return false;
    }

    *isect = SurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv,
        Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this);
    isect->n = isect->shading.n = Normal3f(dp02.cross(dp12).normalized());

    if (mesh->n)
        isect->n = Faceforward(isect->n, isect->shading.n);
    else if (reverseOrientation ^ transformSwapsHandedness)
        isect->n = isect->shading.n = -isect->n;

    Normal3f ns;
    // pre-face normal
    if (mesh->n) ns = (b0 * mesh->n[v[0]] +
        b1 * mesh->n[v[1]] +
        b2 * mesh->n[v[2]]).normalized();
    else
        ns = isect->n;
    
    // shading tangent
    Vector3f ss;
    if (mesh->s) ss = (b0 * mesh->s[v[0]] +
        b1 * mesh->s[v[1]] +
        b2 * mesh->s[v[2]]).normalized();
    else
        ss = (isect->dpdu).normalized();

    // shading bitangent 
    Vector3f ts = ss.cross(ns);
    if (ts.squaredNorm() > 0.f) {
        ts = ts.normalized();
        ss = ts.cross(ns);
    }
    else
        CoordinateSystem((Vector3f)ns, &ss, &ts);


    *tHit = t;
    return true;
}

bool Triangle::IntersectP(const Ray& r, bool testAlphaTexture) const
{

    // TODO:

    return false;
}

void Triangle::GetUVs(Point2f uv[3]) const {
    if (mesh->uv) {
        uv[0] = mesh->uv[v[0]];
        uv[1] = mesh->uv[v[1]];
        uv[2] = mesh->uv[v[2]];
    }
    else {
        uv[0] = Point2f(0, 0);
        uv[1] = Point2f(1, 0);
        uv[2] = Point2f(1, 1);
    }
}

double Triangle::Area() const
{
    const Point3f& p0 = mesh->p[v[0]];
    const Point3f& p1 = mesh->p[v[1]];
    const Point3f& p2 = mesh->p[v[2]];

    return 0.5 * std::sqrt((p1 - p0).cross(p2 - p0).squaredNorm());

}


TriangleMesh::TriangleMesh(const Transform& ObjectToWorld, 
    int nTriangles, 
    const int* vertexIndices,
    int nVertices,
    const Point3f* P, 
    const Vector3f* S,
    const Normal3f* N, 
    const Point2f* UV, 
    const std::shared_ptr<Texture<Float>>& alphaMask) :
    nTriangles(nTriangles), nVertices(nVertices),
    vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles),
    alphaMask(alphaMask)
{
    p.reset(new Point3f[nVertices]);
    for (int i = 0; i < nVertices; ++i)
        p[i] = ObjectToWorld(P[i]);

    if (UV) 
    {
        uv.reset(new Point2f[nVertices]);
        memcpy(uv.get(), UV, nVertices * sizeof(Point2f));
    }
    if (N) 
    {
        n.reset(new Normal3f[nVertices]);
        for (int i = 0; i < nVertices; ++i) 
            n[i] = ObjectToWorld(N[i]);
    }
    if (S) 
    {
        s.reset(new Vector3f[nVertices]);
        for (int i = 0; i < nVertices; ++i)
            s[i] = ObjectToWorld(S[i]);
    }


}
