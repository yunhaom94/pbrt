
#include "core/transform.h"
#include "core/interaction.h"
#include "core/ray.h"
#include "core/bounding_boxes.h"



Transform Inverse(const Transform& t)
{
	return Transform(t.mInv, t.m);
}

Transform Transpose(const Transform& t)
{
	return Transform(t.m.transpose(), t.mInv.transpose());
}

Matrix4x4 Transpose(const Matrix4x4& m) 
{
	Matrix4x4 out;
	out << m.row(0)[0], m.row(1)[0], m.row(2)[0], m.row(3)[0], m.row(0)[1],
		m.row(1)[1], m.row(2)[1], m.row(3)[1], m.row(0)[2], m.row(1)[2],
		m.row(2)[2], m.row(3)[2], m.row(0)[3], m.row(1)[3], m.row(2)[3],
		m.row(3)[3];

	return out;
}

Transform Translate(const Vector3f& delta)
{
	Matrix4x4 m;
	m << 1, 0, 0, delta.x(),
		 0, 1, 0, delta.y(),
		 0, 0, 1, delta.z(),
		 0, 0, 0, 1;

	Matrix4x4 minv;
	minv << 1, 0, 0, -delta.x(),
		    0, 1, 0, -delta.y(),
		    0, 0, 1, -delta.z(),
		    0, 0, 0, 1;

	return Transform(m, minv);
}

Transform Scale(Float x, Float y, Float z)
{
	Matrix4x4 m;
	m << x, 0, 0, 0,
		 0, y, 0, 0,
		 0, 0, z, 0,
		 0, 0, 0, 1;

	Matrix4x4 minv;
    minv << (1.0 / x), 0, 0, 0,
		 0, (1.0 / y), 0, 0,
		 0, 0, (1.0 / z), 0,
		 0, 0, 0, 1;

	return Transform(m, minv);
}

Transform RotateX(Float theta)
{
	Matrix4x4 m;
	m << 1, 0,               0,                0,
		 0, std::cos(theta), -std::sin(theta), 0,
		 0, std::sin(theta), std::cos(theta),  0,
		 0, 0,               0,                1;

	return Transform(m, m.transpose());
}

Transform RotateY(Float theta)
{
	Matrix4x4 m;
	m << std::cos(theta),  0, std::sin(theta), 0,
		 0,                1, 0,               0,
		 -std::sin(theta), 0, std::cos(theta), 0,
		 0,                0, 0,               1;


	return Transform(m, m.transpose());
}

Transform RotateZ(Float theta)
{
	Matrix4x4 m;
	m << std::cos(theta), -std::sin(theta), 0, 0,
		 std::sin(theta), std::cos(theta),  0, 0,
		 0,               0,                1, 0,
		 0,               0,                0, 1;


	return Transform(m, m.transpose());
}

Transform Rotate(Float theta, const Vector3f& axis)
{

	Vector3f a = axis.normalized();
	Float sinTheta = std::sin(theta);
	Float cosTheta = std::cos(theta);
	Matrix4x4 m = Matrix4x4::Zero();

	m.row(0)[0] = a.x() * a.x() + (1 - a.x() * a.x()) * cosTheta;
	m.row(0)[1] = a.x() * a.y() * (1 - cosTheta) - a.z() * sinTheta;
	m.row(0)[2] = a.x() * a.z() * (1 - cosTheta) + a.y() * sinTheta;
	m.row(0)[3] = 0;

	m.row(1)[0] = a.x() * a.y() * (1 - cosTheta) + a.z() * sinTheta;
	m.row(1)[1] = a.y() * a.y() + (1 - a.y() * a.y()) * cosTheta;
	m.row(1)[2] = a.y() * a.z() * (1 - cosTheta) - a.x() * sinTheta;
	m.row(1)[3] = 0;

	m.row(1)[0] = a.x() * a.z() * (1 - cosTheta) - a.y() * sinTheta;
	m.row(1)[1] = a.y() * a.z() * (1 - cosTheta) + a.x() * sinTheta;
	m.row(1)[2] = a.z() * a.z() + (1 - a.z() * a.z()) * cosTheta;
	m.row(1)[3] = 0;

	return Transform();
}

Transform Orthographic(Float zNear, Float zFar)
{
	return Scale(1, 1, 1 / (zFar - zNear)) *
		Translate(Vector3f(0, 0, -zNear));
}

Transform Perspective(Float fov, Float n, Float f)
{

	Matrix4x4 persp;
	persp << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, f / (f - n), -f * n / (f - n),
		0, 0, 1, 0;

	Float invTanAng = 1 / std::tan(Radians(fov) / 2);
	return Scale(invTanAng, invTanAng, 1) * Transform(persp);
}

Transform LookAt(const Vector3f& pos, const Vector3f& look, const Vector3f& up)
{

	Vector3f dir = (look - pos).normalized();
	Vector3f left = (up.normalized().cross(dir)).normalized();
	Vector3f newUp = dir.cross(left);

	Matrix4x4 cameraToWorld;

	cameraToWorld << 
		left.x(), newUp.x(), dir.x(), pos.x(),
		left.y(), newUp.y(), dir.y(), pos.y(),
		left.z(), newUp.z(), dir.z(), pos.z(),
		0,        0,         0,       1;

	return Transform(cameraToWorld.inverse(), cameraToWorld);
}


Ray Transform::operator()(const Ray& r) const
{
	Vector3f oError;
	Point3f o3 = (*this)(r.o);
	Vector3f d3 = (*this)(r.d);

	// TODO: Offset ray origin to edge of error bounds and compute tMax p233

	return Ray(o3, d3, r.tMax, r.time, r.medium);
}

Ray Transform::operator()(const Ray& r, Vector3f* oError, Vector3f* dError) const
{
	Point3f o = (*this)(r.o, oError);
	Vector3f d = (*this)(r.d, dError);

	// TODO: Offset ray origin to edge of error bounds and compute tMax p233

	return Ray(o, d, r.tMax, r.time, r.medium);
}

Ray Transform::operator()(const Ray& r, const Vector3f& oErrorIn, const Vector3f& dErrorIn, Vector3f* oErrorOut, Vector3f* dErrorOut) const
{
	Point3f o = (*this)(r.o, oErrorIn, oErrorOut);
	Vector3f d = (*this)(r.d, dErrorIn, dErrorOut);

	// TODO: Offset ray origin to edge of error bounds and compute tMax p233

	return Ray(o, d, r.tMax, r.time, r.medium);
}

RayDifferential Transform::operator()(const RayDifferential& r) const
{
    Ray tr = (*this)(Ray(r));
    RayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
    ret.hasDifferentials = r.hasDifferentials;
    ret.rxOrigin = (*this)(r.rxOrigin);
    ret.ryOrigin = (*this)(r.ryOrigin);
    ret.rxDirection = (*this)(r.rxDirection);
    ret.ryDirection = (*this)(r.ryDirection);
    return ret;
}


Bounds3f Transform::operator()(const Bounds3f& b) const
{
	const Transform& M = *this;
	Bounds3f ret(M(Point3f(b.pMin.x(), b.pMin.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMin.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMin.x(), b.pMax.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMin.x(), b.pMin.y(), b.pMax.z())));
	ret = Union(ret, M(Point3f(b.pMin.x(), b.pMax.y(), b.pMax.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMax.y(), b.pMin.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMin.y(), b.pMax.z())));
	ret = Union(ret, M(Point3f(b.pMax.x(), b.pMax.y(), b.pMax.z())));
	return ret;


}

Transform Transform::operator*(const Transform& t2) const {
	return Transform(m * t2.m, mInv * t2.mInv);
}

bool Transform::operator==(const Transform& t) const
{
	return t.m == m && t.mInv == mInv;
}

bool Transform::operator!=(const Transform& t) const
{
	return t.m != m || t.mInv != mInv;
}

bool Transform::HasScale() const
{
	Float la2 = (*this)(Vector3f(1, 0, 0)).squaredNorm();
	Float lb2 = (*this)(Vector3f(0, 1, 0)).squaredNorm();
	Float lc2 = (*this)(Vector3f(0, 0, 1)).squaredNorm();
#define NOT_ONE(x) ((x) < 0.999 || (x) > 1.001)
	return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
}

bool Transform::SwapsHandedness() const
{
	Eigen::Matrix3d mtemp;

	mtemp.row(0) = m.row(0).head<3>();
	mtemp.row(1) = m.row(1).head<3>();
	mtemp.row(2) = m.row(2).head<3>();

	return mtemp.determinant() < 0;
}

SurfaceInteraction Transform::operator()(const SurfaceInteraction& si) const
{
	SurfaceInteraction ret;
	// TODO: Transform p and pError in SurfaceInteraction p229

	// TODO: maybe write a clone function
	ret.wo = (*this)(si.wo);
	ret.n = (*this)(si.n);
	ret.time = si.time;
	ret.bsdf = si.bsdf;
	ret.surfaceInteraction = si.surfaceInteraction;
	ret.mediumInterface = si.mediumInterface;
	ret.uv = si.uv;
	ret.dpdu = (*this)(si.dpdu);
	ret.dpdv = (*this)(si.dpdv);
	ret.dndu = (*this)(si.dndu);
	ret.dndv = (*this)(si.dndv);
	ret.shape = si.shape;
	ret.shading.n = (*this)(si.shading.n).normalized();
	ret.shading.dpdu = (*this)(si.shading.dpdu);
	ret.shading.dpdv = (*this)(si.shading.dpdv);
	ret.shading.dndu = (*this)(si.shading.dndu);
	ret.shading.dndv = (*this)(si.shading.dndv);

	return ret;
}

class Interval
{
public:
	// Interval Public Methods
	Interval(Float v) : low(v), high(v) {}
	Interval(Float v0, Float v1)
		: low(std::min(v0, v1)), high(std::max(v0, v1)) {}
	Interval operator+(const Interval& i) const {
		return Interval(low + i.low, high + i.high);
	}
	Interval operator-(const Interval& i) const {
		return Interval(low - i.high, high - i.low);
	}
	Interval operator*(const Interval& i) const {
		return Interval(std::min(std::min(low * i.low, high * i.low),
			std::min(low * i.high, high * i.high)),
			std::max(std::max(low * i.low, high * i.low),
				std::max(low * i.high, high * i.high)));
	}
	Float low, high;
};

inline Interval Sin(const Interval& i)
{

	Float sinLow = std::sin(i.low), sinHigh = std::sin(i.high);
	if (sinLow > sinHigh) std::swap(sinLow, sinHigh);
	if (i.low < Pi / 2 && i.high > Pi / 2) sinHigh = 1.;
	if (i.low < (3.f / 2.f) * Pi && i.high >(3.f / 2.f) * Pi) sinLow = -1.;
	return Interval(sinLow, sinHigh);
}

inline Interval Cos(const Interval& i)
{

	Float cosLow = std::cos(i.low), cosHigh = std::cos(i.high);
	if (cosLow > cosHigh) std::swap(cosLow, cosHigh);
	if (i.low < Pi && i.high > Pi) cosLow = -1.;
	return Interval(cosLow, cosHigh);
}

void IntervalFindZeros(Float c1, Float c2, Float c3, Float c4, Float c5,
	Float theta, Interval tInterval, Float* zeros,
	int* zeroCount, int depth = 8) {
	// Evaluate motion derivative in interval form, return if no zeros
	Interval range = Interval(c1) +
		(Interval(c2) + Interval(c3) * tInterval) *
		Cos(Interval(2 * theta) * tInterval) +
		(Interval(c4) + Interval(c5) * tInterval) *
		Sin(Interval(2 * theta) * tInterval);
	if (range.low > 0. || range.high < 0. || range.low == range.high) return;
	if (depth > 0) {
		// Split _tInterval_ and check both resulting intervals
		Float mid = (tInterval.low + tInterval.high) * 0.5f;
		IntervalFindZeros(c1, c2, c3, c4, c5, theta,
			Interval(tInterval.low, mid), zeros, zeroCount,
			depth - 1);
		IntervalFindZeros(c1, c2, c3, c4, c5, theta,
			Interval(mid, tInterval.high), zeros, zeroCount,
			depth - 1);
	}
	else {
		// Use Newton's method to refine zero
		Float tNewton = (tInterval.low + tInterval.high) * 0.5f;
		for (int i = 0; i < 4; ++i) {
			Float fNewton =
				c1 + (c2 + c3 * tNewton) * std::cos(2.f * theta * tNewton) +
				(c4 + c5 * tNewton) * std::sin(2.f * theta * tNewton);
			Float fPrimeNewton = (c3 + 2 * (c4 + c5 * tNewton) * theta) *
				std::cos(2.f * tNewton * theta) +
				(c5 - 2 * (c2 + c3 * tNewton) * theta) *
				std::sin(2.f * tNewton * theta);
			if (fNewton == 0 || fPrimeNewton == 0) break;
			tNewton = tNewton - fNewton / fPrimeNewton;
		}
		if (tNewton >= tInterval.low - 1e-3f &&
			tNewton < tInterval.high + 1e-3f) {
			zeros[*zeroCount] = tNewton;
			(*zeroCount)++;
		}
	}
}




AnimatedTransform::AnimatedTransform(const Transform* startTransform,
    Float startTime,
    const Transform* endTransform,
    Float endTime)
    : startTransform(startTransform),
    endTransform(endTransform),
    startTime(startTime),
    endTime(endTime),
    actuallyAnimated(*startTransform != *endTransform) {
    if (!actuallyAnimated)
        return;
    Decompose(startTransform->m, &T[0], &R[0], &S[0]);
    Decompose(endTransform->m, &T[1], &R[1], &S[1]);
    // Flip _R[1]_ if needed to select shortest path
    if (Dot(R[0], R[1]) < 0) R[1] = -R[1];
    hasRotation = Dot(R[0], R[1]) < 0.9995f;
    // Compute terms of motion derivative function
    if (hasRotation) {
        Float cosTheta = Dot(R[0], R[1]);
        Float theta = std::acos(Clamp(cosTheta, -1, 1));
        Quaternion qperp = Normalize(R[1] - R[0] * cosTheta);

        Float t0x = T[0].x();
        Float t0y = T[0].y();
        Float t0z = T[0].z();
        Float t1x = T[1].x();
        Float t1y = T[1].y();
        Float t1z = T[1].z();
        Float q0x = R[0].v.x();
        Float q0y = R[0].v.y();
        Float q0z = R[0].v.z();
        Float q0w = R[0].w;
        Float qperpx = qperp.v.x();
        Float qperpy = qperp.v.y();
        Float qperpz = qperp.v.z();
        Float qperpw = qperp.w;
        Float s000 = S[0].row(0)[0];
        Float s001 = S[0].row(0)[1];
        Float s002 = S[0].row(0)[2];
        Float s010 = S[0].row(1)[0];
        Float s011 = S[0].row(1)[1];
        Float s012 = S[0].row(1)[2];
        Float s020 = S[0].row(2)[0];
        Float s021 = S[0].row(2)[1];
        Float s022 = S[0].row(2)[2];
        Float s100 = S[1].row(0)[0];
        Float s101 = S[1].row(0)[1];
        Float s102 = S[1].row(0)[2];
        Float s110 = S[1].row(1)[0];
        Float s111 = S[1].row(1)[1];
        Float s112 = S[1].row(1)[2];
        Float s120 = S[1].row(2)[0];
        Float s121 = S[1].row(2)[1];
        Float s122 = S[1].row(2)[2];

        // Compute terms of motion derivative function
        c1[0] = DerivativeTerm(
            -t0x + t1x,
            (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
            s000 +
            q0w * q0z * s010 - qperpx * qperpy * s010 +
            qperpw * qperpz * s010 - q0w * q0y * s020 -
            qperpw * qperpy * s020 - qperpx * qperpz * s020 + s100 -
            q0y * q0y * s100 - q0z * q0z * s100 - qperpy * qperpy * s100 -
            qperpz * qperpz * s100 - q0w * q0z * s110 +
            qperpx * qperpy * s110 - qperpw * qperpz * s110 +
            q0w * q0y * s120 + qperpw * qperpy * s120 +
            qperpx * qperpz * s120 +
            q0x * (-(q0y * s010) - q0z * s020 + q0y * s110 + q0z * s120),
            (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
            s001 +
            q0w * q0z * s011 - qperpx * qperpy * s011 +
            qperpw * qperpz * s011 - q0w * q0y * s021 -
            qperpw * qperpy * s021 - qperpx * qperpz * s021 + s101 -
            q0y * q0y * s101 - q0z * q0z * s101 - qperpy * qperpy * s101 -
            qperpz * qperpz * s101 - q0w * q0z * s111 +
            qperpx * qperpy * s111 - qperpw * qperpz * s111 +
            q0w * q0y * s121 + qperpw * qperpy * s121 +
            qperpx * qperpz * s121 +
            q0x * (-(q0y * s011) - q0z * s021 + q0y * s111 + q0z * s121),
            (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
            s002 +
            q0w * q0z * s012 - qperpx * qperpy * s012 +
            qperpw * qperpz * s012 - q0w * q0y * s022 -
            qperpw * qperpy * s022 - qperpx * qperpz * s022 + s102 -
            q0y * q0y * s102 - q0z * q0z * s102 - qperpy * qperpy * s102 -
            qperpz * qperpz * s102 - q0w * q0z * s112 +
            qperpx * qperpy * s112 - qperpw * qperpz * s112 +
            q0w * q0y * s122 + qperpw * qperpy * s122 +
            qperpx * qperpz * s122 +
            q0x * (-(q0y * s012) - q0z * s022 + q0y * s112 + q0z * s122));

        c2[0] = DerivativeTerm(
            0.,
            -(qperpy * qperpy * s000) - qperpz * qperpz * s000 +
            qperpx * qperpy * s010 - qperpw * qperpz * s010 +
            qperpw * qperpy * s020 + qperpx * qperpz * s020 +
            q0y * q0y * (s000 - s100) + q0z * q0z * (s000 - s100) +
            qperpy * qperpy * s100 + qperpz * qperpz * s100 -
            qperpx * qperpy * s110 + qperpw * qperpz * s110 -
            qperpw * qperpy * s120 - qperpx * qperpz * s120 +
            2 * q0x * qperpy * s010 * theta -
            2 * q0w * qperpz * s010 * theta +
            2 * q0w * qperpy * s020 * theta +
            2 * q0x * qperpz * s020 * theta +
            q0y *
            (q0x * (-s010 + s110) + q0w * (-s020 + s120) +
                2 * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020) *
                theta) +
            q0z * (q0w * (s010 - s110) + q0x * (-s020 + s120) -
                2 * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020) *
                theta),
            -(qperpy * qperpy * s001) - qperpz * qperpz * s001 +
            qperpx * qperpy * s011 - qperpw * qperpz * s011 +
            qperpw * qperpy * s021 + qperpx * qperpz * s021 +
            q0y * q0y * (s001 - s101) + q0z * q0z * (s001 - s101) +
            qperpy * qperpy * s101 + qperpz * qperpz * s101 -
            qperpx * qperpy * s111 + qperpw * qperpz * s111 -
            qperpw * qperpy * s121 - qperpx * qperpz * s121 +
            2 * q0x * qperpy * s011 * theta -
            2 * q0w * qperpz * s011 * theta +
            2 * q0w * qperpy * s021 * theta +
            2 * q0x * qperpz * s021 * theta +
            q0y *
            (q0x * (-s011 + s111) + q0w * (-s021 + s121) +
                2 * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021) *
                theta) +
            q0z * (q0w * (s011 - s111) + q0x * (-s021 + s121) -
                2 * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021) *
                theta),
            -(qperpy * qperpy * s002) - qperpz * qperpz * s002 +
            qperpx * qperpy * s012 - qperpw * qperpz * s012 +
            qperpw * qperpy * s022 + qperpx * qperpz * s022 +
            q0y * q0y * (s002 - s102) + q0z * q0z * (s002 - s102) +
            qperpy * qperpy * s102 + qperpz * qperpz * s102 -
            qperpx * qperpy * s112 + qperpw * qperpz * s112 -
            qperpw * qperpy * s122 - qperpx * qperpz * s122 +
            2 * q0x * qperpy * s012 * theta -
            2 * q0w * qperpz * s012 * theta +
            2 * q0w * qperpy * s022 * theta +
            2 * q0x * qperpz * s022 * theta +
            q0y *
            (q0x * (-s012 + s112) + q0w * (-s022 + s122) +
                2 * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022) *
                theta) +
            q0z * (q0w * (s012 - s112) + q0x * (-s022 + s122) -
                2 * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022) *
                theta));

        c3[0] = DerivativeTerm(
            0.,
            -2 * (q0x * qperpy * s010 - q0w * qperpz * s010 +
                q0w * qperpy * s020 + q0x * qperpz * s020 -
                q0x * qperpy * s110 + q0w * qperpz * s110 -
                q0w * qperpy * s120 - q0x * qperpz * s120 +
                q0y * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020 +
                    2 * qperpy * s100 - qperpx * s110 - qperpw * s120) +
                q0z * (-2 * qperpz * s000 - qperpw * s010 + qperpx * s020 +
                    2 * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
            theta,
            -2 * (q0x * qperpy * s011 - q0w * qperpz * s011 +
                q0w * qperpy * s021 + q0x * qperpz * s021 -
                q0x * qperpy * s111 + q0w * qperpz * s111 -
                q0w * qperpy * s121 - q0x * qperpz * s121 +
                q0y * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021 +
                    2 * qperpy * s101 - qperpx * s111 - qperpw * s121) +
                q0z * (-2 * qperpz * s001 - qperpw * s011 + qperpx * s021 +
                    2 * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
            theta,
            -2 * (q0x * qperpy * s012 - q0w * qperpz * s012 +
                q0w * qperpy * s022 + q0x * qperpz * s022 -
                q0x * qperpy * s112 + q0w * qperpz * s112 -
                q0w * qperpy * s122 - q0x * qperpz * s122 +
                q0y * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022 +
                    2 * qperpy * s102 - qperpx * s112 - qperpw * s122) +
                q0z * (-2 * qperpz * s002 - qperpw * s012 + qperpx * s022 +
                    2 * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
            theta);

        c4[0] = DerivativeTerm(
            0.,
            -(q0x * qperpy * s010) + q0w * qperpz * s010 - q0w * qperpy * s020 -
            q0x * qperpz * s020 + q0x * qperpy * s110 -
            q0w * qperpz * s110 + q0w * qperpy * s120 +
            q0x * qperpz * s120 + 2 * q0y * q0y * s000 * theta +
            2 * q0z * q0z * s000 * theta -
            2 * qperpy * qperpy * s000 * theta -
            2 * qperpz * qperpz * s000 * theta +
            2 * qperpx * qperpy * s010 * theta -
            2 * qperpw * qperpz * s010 * theta +
            2 * qperpw * qperpy * s020 * theta +
            2 * qperpx * qperpz * s020 * theta +
            q0y * (-(qperpx * s010) - qperpw * s020 +
                2 * qperpy * (s000 - s100) + qperpx * s110 +
                qperpw * s120 - 2 * q0x * s010 * theta -
                2 * q0w * s020 * theta) +
            q0z * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020 -
                2 * qperpz * s100 - qperpw * s110 + qperpx * s120 +
                2 * q0w * s010 * theta - 2 * q0x * s020 * theta),
            -(q0x * qperpy * s011) + q0w * qperpz * s011 - q0w * qperpy * s021 -
            q0x * qperpz * s021 + q0x * qperpy * s111 -
            q0w * qperpz * s111 + q0w * qperpy * s121 +
            q0x * qperpz * s121 + 2 * q0y * q0y * s001 * theta +
            2 * q0z * q0z * s001 * theta -
            2 * qperpy * qperpy * s001 * theta -
            2 * qperpz * qperpz * s001 * theta +
            2 * qperpx * qperpy * s011 * theta -
            2 * qperpw * qperpz * s011 * theta +
            2 * qperpw * qperpy * s021 * theta +
            2 * qperpx * qperpz * s021 * theta +
            q0y * (-(qperpx * s011) - qperpw * s021 +
                2 * qperpy * (s001 - s101) + qperpx * s111 +
                qperpw * s121 - 2 * q0x * s011 * theta -
                2 * q0w * s021 * theta) +
            q0z * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021 -
                2 * qperpz * s101 - qperpw * s111 + qperpx * s121 +
                2 * q0w * s011 * theta - 2 * q0x * s021 * theta),
            -(q0x * qperpy * s012) + q0w * qperpz * s012 - q0w * qperpy * s022 -
            q0x * qperpz * s022 + q0x * qperpy * s112 -
            q0w * qperpz * s112 + q0w * qperpy * s122 +
            q0x * qperpz * s122 + 2 * q0y * q0y * s002 * theta +
            2 * q0z * q0z * s002 * theta -
            2 * qperpy * qperpy * s002 * theta -
            2 * qperpz * qperpz * s002 * theta +
            2 * qperpx * qperpy * s012 * theta -
            2 * qperpw * qperpz * s012 * theta +
            2 * qperpw * qperpy * s022 * theta +
            2 * qperpx * qperpz * s022 * theta +
            q0y * (-(qperpx * s012) - qperpw * s022 +
                2 * qperpy * (s002 - s102) + qperpx * s112 +
                qperpw * s122 - 2 * q0x * s012 * theta -
                2 * q0w * s022 * theta) +
            q0z * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022 -
                2 * qperpz * s102 - qperpw * s112 + qperpx * s122 +
                2 * q0w * s012 * theta - 2 * q0x * s022 * theta));

        c5[0] = DerivativeTerm(
            0.,
            2 * (qperpy * qperpy * s000 + qperpz * qperpz * s000 -
                qperpx * qperpy * s010 + qperpw * qperpz * s010 -
                qperpw * qperpy * s020 - qperpx * qperpz * s020 -
                qperpy * qperpy * s100 - qperpz * qperpz * s100 +
                q0y * q0y * (-s000 + s100) + q0z * q0z * (-s000 + s100) +
                qperpx * qperpy * s110 - qperpw * qperpz * s110 +
                q0y * (q0x * (s010 - s110) + q0w * (s020 - s120)) +
                qperpw * qperpy * s120 + qperpx * qperpz * s120 +
                q0z * (-(q0w * s010) + q0x * s020 + q0w * s110 - q0x * s120)) *
            theta,
            2 * (qperpy * qperpy * s001 + qperpz * qperpz * s001 -
                qperpx * qperpy * s011 + qperpw * qperpz * s011 -
                qperpw * qperpy * s021 - qperpx * qperpz * s021 -
                qperpy * qperpy * s101 - qperpz * qperpz * s101 +
                q0y * q0y * (-s001 + s101) + q0z * q0z * (-s001 + s101) +
                qperpx * qperpy * s111 - qperpw * qperpz * s111 +
                q0y * (q0x * (s011 - s111) + q0w * (s021 - s121)) +
                qperpw * qperpy * s121 + qperpx * qperpz * s121 +
                q0z * (-(q0w * s011) + q0x * s021 + q0w * s111 - q0x * s121)) *
            theta,
            2 * (qperpy * qperpy * s002 + qperpz * qperpz * s002 -
                qperpx * qperpy * s012 + qperpw * qperpz * s012 -
                qperpw * qperpy * s022 - qperpx * qperpz * s022 -
                qperpy * qperpy * s102 - qperpz * qperpz * s102 +
                q0y * q0y * (-s002 + s102) + q0z * q0z * (-s002 + s102) +
                qperpx * qperpy * s112 - qperpw * qperpz * s112 +
                q0y * (q0x * (s012 - s112) + q0w * (s022 - s122)) +
                qperpw * qperpy * s122 + qperpx * qperpz * s122 +
                q0z * (-(q0w * s012) + q0x * s022 + q0w * s112 - q0x * s122)) *
            theta);

        c1[1] = DerivativeTerm(
            -t0y + t1y,
            -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 +
            q0z * q0z * s010 + qperpx * qperpx * s010 +
            qperpz * qperpz * s010 - q0y * q0z * s020 +
            qperpw * qperpx * s020 - qperpy * qperpz * s020 +
            qperpx * qperpy * s100 + qperpw * qperpz * s100 +
            q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) + s110 -
            q0z * q0z * s110 - qperpx * qperpx * s110 -
            qperpz * qperpz * s110 +
            q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
            q0y * q0z * s120 - qperpw * qperpx * s120 +
            qperpy * qperpz * s120,
            -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 +
            q0z * q0z * s011 + qperpx * qperpx * s011 +
            qperpz * qperpz * s011 - q0y * q0z * s021 +
            qperpw * qperpx * s021 - qperpy * qperpz * s021 +
            qperpx * qperpy * s101 + qperpw * qperpz * s101 +
            q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) + s111 -
            q0z * q0z * s111 - qperpx * qperpx * s111 -
            qperpz * qperpz * s111 +
            q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
            q0y * q0z * s121 - qperpw * qperpx * s121 +
            qperpy * qperpz * s121,
            -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 +
            q0z * q0z * s012 + qperpx * qperpx * s012 +
            qperpz * qperpz * s012 - q0y * q0z * s022 +
            qperpw * qperpx * s022 - qperpy * qperpz * s022 +
            qperpx * qperpy * s102 + qperpw * qperpz * s102 +
            q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) + s112 -
            q0z * q0z * s112 - qperpx * qperpx * s112 -
            qperpz * qperpz * s112 +
            q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
            q0y * q0z * s122 - qperpw * qperpx * s122 +
            qperpy * qperpz * s122);

        c2[1] = DerivativeTerm(
            0.,
            qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010 -
            qperpx * qperpx * s010 - qperpz * qperpz * s010 -
            q0y * q0z * s020 - qperpw * qperpx * s020 +
            qperpy * qperpz * s020 - qperpx * qperpy * s100 -
            qperpw * qperpz * s100 + q0x * q0x * (s010 - s110) -
            q0z * q0z * s110 + qperpx * qperpx * s110 +
            qperpz * qperpz * s110 + q0y * q0z * s120 +
            qperpw * qperpx * s120 - qperpy * qperpz * s120 +
            2 * q0z * qperpw * s000 * theta +
            2 * q0y * qperpx * s000 * theta -
            4 * q0z * qperpz * s010 * theta +
            2 * q0z * qperpy * s020 * theta +
            2 * q0y * qperpz * s020 * theta +
            q0x * (q0w * s020 + q0y * (-s000 + s100) - q0w * s120 +
                2 * qperpy * s000 * theta - 4 * qperpx * s010 * theta -
                2 * qperpw * s020 * theta) +
            q0w * (-(q0z * s000) + q0z * s100 + 2 * qperpz * s000 * theta -
                2 * qperpx * s020 * theta),
            qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011 -
            qperpx * qperpx * s011 - qperpz * qperpz * s011 -
            q0y * q0z * s021 - qperpw * qperpx * s021 +
            qperpy * qperpz * s021 - qperpx * qperpy * s101 -
            qperpw * qperpz * s101 + q0x * q0x * (s011 - s111) -
            q0z * q0z * s111 + qperpx * qperpx * s111 +
            qperpz * qperpz * s111 + q0y * q0z * s121 +
            qperpw * qperpx * s121 - qperpy * qperpz * s121 +
            2 * q0z * qperpw * s001 * theta +
            2 * q0y * qperpx * s001 * theta -
            4 * q0z * qperpz * s011 * theta +
            2 * q0z * qperpy * s021 * theta +
            2 * q0y * qperpz * s021 * theta +
            q0x * (q0w * s021 + q0y * (-s001 + s101) - q0w * s121 +
                2 * qperpy * s001 * theta - 4 * qperpx * s011 * theta -
                2 * qperpw * s021 * theta) +
            q0w * (-(q0z * s001) + q0z * s101 + 2 * qperpz * s001 * theta -
                2 * qperpx * s021 * theta),
            qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012 -
            qperpx * qperpx * s012 - qperpz * qperpz * s012 -
            q0y * q0z * s022 - qperpw * qperpx * s022 +
            qperpy * qperpz * s022 - qperpx * qperpy * s102 -
            qperpw * qperpz * s102 + q0x * q0x * (s012 - s112) -
            q0z * q0z * s112 + qperpx * qperpx * s112 +
            qperpz * qperpz * s112 + q0y * q0z * s122 +
            qperpw * qperpx * s122 - qperpy * qperpz * s122 +
            2 * q0z * qperpw * s002 * theta +
            2 * q0y * qperpx * s002 * theta -
            4 * q0z * qperpz * s012 * theta +
            2 * q0z * qperpy * s022 * theta +
            2 * q0y * qperpz * s022 * theta +
            q0x * (q0w * s022 + q0y * (-s002 + s102) - q0w * s122 +
                2 * qperpy * s002 * theta - 4 * qperpx * s012 * theta -
                2 * qperpw * s022 * theta) +
            q0w * (-(q0z * s002) + q0z * s102 + 2 * qperpz * s002 * theta -
                2 * qperpx * s022 * theta));

        c3[1] = DerivativeTerm(
            0., 2 * (-(q0x * qperpy * s000) - q0w * qperpz * s000 +
                2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
                q0w * qperpx * s020 + q0x * qperpy * s100 +
                q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
                q0x * qperpw * s120 - q0w * qperpx * s120 +
                q0z * (2 * qperpz * s010 - qperpy * s020 +
                    qperpw * (-s000 + s100) - 2 * qperpz * s110 +
                    qperpy * s120) +
                q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                    qperpz * s120)) *
            theta,
            2 * (-(q0x * qperpy * s001) - q0w * qperpz * s001 +
                2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
                q0w * qperpx * s021 + q0x * qperpy * s101 +
                q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
                q0x * qperpw * s121 - q0w * qperpx * s121 +
                q0z * (2 * qperpz * s011 - qperpy * s021 +
                    qperpw * (-s001 + s101) - 2 * qperpz * s111 +
                    qperpy * s121) +
                q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                    qperpz * s121)) *
            theta,
            2 * (-(q0x * qperpy * s002) - q0w * qperpz * s002 +
                2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
                q0w * qperpx * s022 + q0x * qperpy * s102 +
                q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
                q0x * qperpw * s122 - q0w * qperpx * s122 +
                q0z * (2 * qperpz * s012 - qperpy * s022 +
                    qperpw * (-s002 + s102) - 2 * qperpz * s112 +
                    qperpy * s122) +
                q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                    qperpz * s122)) *
            theta);

        c4[1] = DerivativeTerm(
            0.,
            -(q0x * qperpy * s000) - q0w * qperpz * s000 +
            2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
            q0w * qperpx * s020 + q0x * qperpy * s100 +
            q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
            q0x * qperpw * s120 - q0w * qperpx * s120 +
            2 * qperpx * qperpy * s000 * theta +
            2 * qperpw * qperpz * s000 * theta +
            2 * q0x * q0x * s010 * theta + 2 * q0z * q0z * s010 * theta -
            2 * qperpx * qperpx * s010 * theta -
            2 * qperpz * qperpz * s010 * theta +
            2 * q0w * q0x * s020 * theta -
            2 * qperpw * qperpx * s020 * theta +
            2 * qperpy * qperpz * s020 * theta +
            q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                qperpz * s120 - 2 * q0x * s000 * theta) +
            q0z * (2 * qperpz * s010 - qperpy * s020 +
                qperpw * (-s000 + s100) - 2 * qperpz * s110 +
                qperpy * s120 - 2 * q0w * s000 * theta -
                2 * q0y * s020 * theta),
            -(q0x * qperpy * s001) - q0w * qperpz * s001 +
            2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
            q0w * qperpx * s021 + q0x * qperpy * s101 +
            q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
            q0x * qperpw * s121 - q0w * qperpx * s121 +
            2 * qperpx * qperpy * s001 * theta +
            2 * qperpw * qperpz * s001 * theta +
            2 * q0x * q0x * s011 * theta + 2 * q0z * q0z * s011 * theta -
            2 * qperpx * qperpx * s011 * theta -
            2 * qperpz * qperpz * s011 * theta +
            2 * q0w * q0x * s021 * theta -
            2 * qperpw * qperpx * s021 * theta +
            2 * qperpy * qperpz * s021 * theta +
            q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                qperpz * s121 - 2 * q0x * s001 * theta) +
            q0z * (2 * qperpz * s011 - qperpy * s021 +
                qperpw * (-s001 + s101) - 2 * qperpz * s111 +
                qperpy * s121 - 2 * q0w * s001 * theta -
                2 * q0y * s021 * theta),
            -(q0x * qperpy * s002) - q0w * qperpz * s002 +
            2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
            q0w * qperpx * s022 + q0x * qperpy * s102 +
            q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
            q0x * qperpw * s122 - q0w * qperpx * s122 +
            2 * qperpx * qperpy * s002 * theta +
            2 * qperpw * qperpz * s002 * theta +
            2 * q0x * q0x * s012 * theta + 2 * q0z * q0z * s012 * theta -
            2 * qperpx * qperpx * s012 * theta -
            2 * qperpz * qperpz * s012 * theta +
            2 * q0w * q0x * s022 * theta -
            2 * qperpw * qperpx * s022 * theta +
            2 * qperpy * qperpz * s022 * theta +
            q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                qperpz * s122 - 2 * q0x * s002 * theta) +
            q0z * (2 * qperpz * s012 - qperpy * s022 +
                qperpw * (-s002 + s102) - 2 * qperpz * s112 +
                qperpy * s122 - 2 * q0w * s002 * theta -
                2 * q0y * s022 * theta));

        c5[1] = DerivativeTerm(
            0., -2 * (qperpx * qperpy * s000 + qperpw * qperpz * s000 +
                q0z * q0z * s010 - qperpx * qperpx * s010 -
                qperpz * qperpz * s010 - q0y * q0z * s020 -
                qperpw * qperpx * s020 + qperpy * qperpz * s020 -
                qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) -
                q0z * q0z * s110 + qperpx * qperpx * s110 +
                qperpz * qperpz * s110 +
                q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
                q0y * q0z * s120 + qperpw * qperpx * s120 -
                qperpy * qperpz * s120) *
            theta,
            -2 * (qperpx * qperpy * s001 + qperpw * qperpz * s001 +
                q0z * q0z * s011 - qperpx * qperpx * s011 -
                qperpz * qperpz * s011 - q0y * q0z * s021 -
                qperpw * qperpx * s021 + qperpy * qperpz * s021 -
                qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) -
                q0z * q0z * s111 + qperpx * qperpx * s111 +
                qperpz * qperpz * s111 +
                q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
                q0y * q0z * s121 + qperpw * qperpx * s121 -
                qperpy * qperpz * s121) *
            theta,
            -2 * (qperpx * qperpy * s002 + qperpw * qperpz * s002 +
                q0z * q0z * s012 - qperpx * qperpx * s012 -
                qperpz * qperpz * s012 - q0y * q0z * s022 -
                qperpw * qperpx * s022 + qperpy * qperpz * s022 -
                qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) -
                q0z * q0z * s112 + qperpx * qperpx * s112 +
                qperpz * qperpz * s112 +
                q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
                q0y * q0z * s122 + qperpw * qperpx * s122 -
                qperpy * qperpz * s122) *
            theta);

        c1[2] = DerivativeTerm(
            -t0z + t1z, (qperpw * qperpy * s000 - qperpx * qperpz * s000 -
                q0y * q0z * s010 - qperpw * qperpx * s010 -
                qperpy * qperpz * s010 - s020 + q0y * q0y * s020 +
                qperpx * qperpx * s020 + qperpy * qperpy * s020 -
                qperpw * qperpy * s100 + qperpx * qperpz * s100 +
                q0x * q0z * (-s000 + s100) + q0y * q0z * s110 +
                qperpw * qperpx * s110 + qperpy * qperpz * s110 +
                q0w * (q0y * (s000 - s100) + q0x * (-s010 + s110)) +
                q0x * q0x * (s020 - s120) + s120 - q0y * q0y * s120 -
                qperpx * qperpx * s120 - qperpy * qperpy * s120),
            (qperpw * qperpy * s001 - qperpx * qperpz * s001 -
                q0y * q0z * s011 - qperpw * qperpx * s011 -
                qperpy * qperpz * s011 - s021 + q0y * q0y * s021 +
                qperpx * qperpx * s021 + qperpy * qperpy * s021 -
                qperpw * qperpy * s101 + qperpx * qperpz * s101 +
                q0x * q0z * (-s001 + s101) + q0y * q0z * s111 +
                qperpw * qperpx * s111 + qperpy * qperpz * s111 +
                q0w * (q0y * (s001 - s101) + q0x * (-s011 + s111)) +
                q0x * q0x * (s021 - s121) + s121 - q0y * q0y * s121 -
                qperpx * qperpx * s121 - qperpy * qperpy * s121),
            (qperpw * qperpy * s002 - qperpx * qperpz * s002 -
                q0y * q0z * s012 - qperpw * qperpx * s012 -
                qperpy * qperpz * s012 - s022 + q0y * q0y * s022 +
                qperpx * qperpx * s022 + qperpy * qperpy * s022 -
                qperpw * qperpy * s102 + qperpx * qperpz * s102 +
                q0x * q0z * (-s002 + s102) + q0y * q0z * s112 +
                qperpw * qperpx * s112 + qperpy * qperpz * s112 +
                q0w * (q0y * (s002 - s102) + q0x * (-s012 + s112)) +
                q0x * q0x * (s022 - s122) + s122 - q0y * q0y * s122 -
                qperpx * qperpx * s122 - qperpy * qperpy * s122));

        c2[2] = DerivativeTerm(
            0.,
            (q0w * q0y * s000 - q0x * q0z * s000 - qperpw * qperpy * s000 +
                qperpx * qperpz * s000 - q0w * q0x * s010 - q0y * q0z * s010 +
                qperpw * qperpx * s010 + qperpy * qperpz * s010 +
                q0x * q0x * s020 + q0y * q0y * s020 - qperpx * qperpx * s020 -
                qperpy * qperpy * s020 - q0w * q0y * s100 + q0x * q0z * s100 +
                qperpw * qperpy * s100 - qperpx * qperpz * s100 +
                q0w * q0x * s110 + q0y * q0z * s110 - qperpw * qperpx * s110 -
                qperpy * qperpz * s110 - q0x * q0x * s120 - q0y * q0y * s120 +
                qperpx * qperpx * s120 + qperpy * qperpy * s120 -
                2 * q0y * qperpw * s000 * theta + 2 * q0z * qperpx * s000 * theta -
                2 * q0w * qperpy * s000 * theta + 2 * q0x * qperpz * s000 * theta +
                2 * q0x * qperpw * s010 * theta + 2 * q0w * qperpx * s010 * theta +
                2 * q0z * qperpy * s010 * theta + 2 * q0y * qperpz * s010 * theta -
                4 * q0x * qperpx * s020 * theta - 4 * q0y * qperpy * s020 * theta),
            (q0w * q0y * s001 - q0x * q0z * s001 - qperpw * qperpy * s001 +
                qperpx * qperpz * s001 - q0w * q0x * s011 - q0y * q0z * s011 +
                qperpw * qperpx * s011 + qperpy * qperpz * s011 +
                q0x * q0x * s021 + q0y * q0y * s021 - qperpx * qperpx * s021 -
                qperpy * qperpy * s021 - q0w * q0y * s101 + q0x * q0z * s101 +
                qperpw * qperpy * s101 - qperpx * qperpz * s101 +
                q0w * q0x * s111 + q0y * q0z * s111 - qperpw * qperpx * s111 -
                qperpy * qperpz * s111 - q0x * q0x * s121 - q0y * q0y * s121 +
                qperpx * qperpx * s121 + qperpy * qperpy * s121 -
                2 * q0y * qperpw * s001 * theta + 2 * q0z * qperpx * s001 * theta -
                2 * q0w * qperpy * s001 * theta + 2 * q0x * qperpz * s001 * theta +
                2 * q0x * qperpw * s011 * theta + 2 * q0w * qperpx * s011 * theta +
                2 * q0z * qperpy * s011 * theta + 2 * q0y * qperpz * s011 * theta -
                4 * q0x * qperpx * s021 * theta - 4 * q0y * qperpy * s021 * theta),
            (q0w * q0y * s002 - q0x * q0z * s002 - qperpw * qperpy * s002 +
                qperpx * qperpz * s002 - q0w * q0x * s012 - q0y * q0z * s012 +
                qperpw * qperpx * s012 + qperpy * qperpz * s012 +
                q0x * q0x * s022 + q0y * q0y * s022 - qperpx * qperpx * s022 -
                qperpy * qperpy * s022 - q0w * q0y * s102 + q0x * q0z * s102 +
                qperpw * qperpy * s102 - qperpx * qperpz * s102 +
                q0w * q0x * s112 + q0y * q0z * s112 - qperpw * qperpx * s112 -
                qperpy * qperpz * s112 - q0x * q0x * s122 - q0y * q0y * s122 +
                qperpx * qperpx * s122 + qperpy * qperpy * s122 -
                2 * q0y * qperpw * s002 * theta + 2 * q0z * qperpx * s002 * theta -
                2 * q0w * qperpy * s002 * theta + 2 * q0x * qperpz * s002 * theta +
                2 * q0x * qperpw * s012 * theta + 2 * q0w * qperpx * s012 * theta +
                2 * q0z * qperpy * s012 * theta + 2 * q0y * qperpz * s012 * theta -
                4 * q0x * qperpx * s022 * theta -
                4 * q0y * qperpy * s022 * theta));

        c3[2] = DerivativeTerm(
            0., -2 * (-(q0w * qperpy * s000) + q0x * qperpz * s000 +
                q0x * qperpw * s010 + q0w * qperpx * s010 -
                2 * q0x * qperpx * s020 + q0w * qperpy * s100 -
                q0x * qperpz * s100 - q0x * qperpw * s110 -
                q0w * qperpx * s110 +
                q0z * (qperpx * s000 + qperpy * s010 - qperpx * s100 -
                    qperpy * s110) +
                2 * q0x * qperpx * s120 +
                q0y * (qperpz * s010 - 2 * qperpy * s020 +
                    qperpw * (-s000 + s100) - qperpz * s110 +
                    2 * qperpy * s120)) *
            theta,
            -2 * (-(q0w * qperpy * s001) + q0x * qperpz * s001 +
                q0x * qperpw * s011 + q0w * qperpx * s011 -
                2 * q0x * qperpx * s021 + q0w * qperpy * s101 -
                q0x * qperpz * s101 - q0x * qperpw * s111 -
                q0w * qperpx * s111 +
                q0z * (qperpx * s001 + qperpy * s011 - qperpx * s101 -
                    qperpy * s111) +
                2 * q0x * qperpx * s121 +
                q0y * (qperpz * s011 - 2 * qperpy * s021 +
                    qperpw * (-s001 + s101) - qperpz * s111 +
                    2 * qperpy * s121)) *
            theta,
            -2 * (-(q0w * qperpy * s002) + q0x * qperpz * s002 +
                q0x * qperpw * s012 + q0w * qperpx * s012 -
                2 * q0x * qperpx * s022 + q0w * qperpy * s102 -
                q0x * qperpz * s102 - q0x * qperpw * s112 -
                q0w * qperpx * s112 +
                q0z * (qperpx * s002 + qperpy * s012 - qperpx * s102 -
                    qperpy * s112) +
                2 * q0x * qperpx * s122 +
                q0y * (qperpz * s012 - 2 * qperpy * s022 +
                    qperpw * (-s002 + s102) - qperpz * s112 +
                    2 * qperpy * s122)) *
            theta);

        c4[2] = DerivativeTerm(
            0.,
            q0w * qperpy * s000 - q0x * qperpz * s000 - q0x * qperpw * s010 -
            q0w * qperpx * s010 + 2 * q0x * qperpx * s020 -
            q0w * qperpy * s100 + q0x * qperpz * s100 +
            q0x * qperpw * s110 + q0w * qperpx * s110 -
            2 * q0x * qperpx * s120 - 2 * qperpw * qperpy * s000 * theta +
            2 * qperpx * qperpz * s000 * theta -
            2 * q0w * q0x * s010 * theta +
            2 * qperpw * qperpx * s010 * theta +
            2 * qperpy * qperpz * s010 * theta +
            2 * q0x * q0x * s020 * theta + 2 * q0y * q0y * s020 * theta -
            2 * qperpx * qperpx * s020 * theta -
            2 * qperpy * qperpy * s020 * theta +
            q0z * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 +
                qperpy * s110 - 2 * q0x * s000 * theta) +
            q0y * (-(qperpz * s010) + 2 * qperpy * s020 +
                qperpw * (s000 - s100) + qperpz * s110 -
                2 * qperpy * s120 + 2 * q0w * s000 * theta -
                2 * q0z * s010 * theta),
            q0w * qperpy * s001 - q0x * qperpz * s001 - q0x * qperpw * s011 -
            q0w * qperpx * s011 + 2 * q0x * qperpx * s021 -
            q0w * qperpy * s101 + q0x * qperpz * s101 +
            q0x * qperpw * s111 + q0w * qperpx * s111 -
            2 * q0x * qperpx * s121 - 2 * qperpw * qperpy * s001 * theta +
            2 * qperpx * qperpz * s001 * theta -
            2 * q0w * q0x * s011 * theta +
            2 * qperpw * qperpx * s011 * theta +
            2 * qperpy * qperpz * s011 * theta +
            2 * q0x * q0x * s021 * theta + 2 * q0y * q0y * s021 * theta -
            2 * qperpx * qperpx * s021 * theta -
            2 * qperpy * qperpy * s021 * theta +
            q0z * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 +
                qperpy * s111 - 2 * q0x * s001 * theta) +
            q0y * (-(qperpz * s011) + 2 * qperpy * s021 +
                qperpw * (s001 - s101) + qperpz * s111 -
                2 * qperpy * s121 + 2 * q0w * s001 * theta -
                2 * q0z * s011 * theta),
            q0w * qperpy * s002 - q0x * qperpz * s002 - q0x * qperpw * s012 -
            q0w * qperpx * s012 + 2 * q0x * qperpx * s022 -
            q0w * qperpy * s102 + q0x * qperpz * s102 +
            q0x * qperpw * s112 + q0w * qperpx * s112 -
            2 * q0x * qperpx * s122 - 2 * qperpw * qperpy * s002 * theta +
            2 * qperpx * qperpz * s002 * theta -
            2 * q0w * q0x * s012 * theta +
            2 * qperpw * qperpx * s012 * theta +
            2 * qperpy * qperpz * s012 * theta +
            2 * q0x * q0x * s022 * theta + 2 * q0y * q0y * s022 * theta -
            2 * qperpx * qperpx * s022 * theta -
            2 * qperpy * qperpy * s022 * theta +
            q0z * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 +
                qperpy * s112 - 2 * q0x * s002 * theta) +
            q0y * (-(qperpz * s012) + 2 * qperpy * s022 +
                qperpw * (s002 - s102) + qperpz * s112 -
                2 * qperpy * s122 + 2 * q0w * s002 * theta -
                2 * q0z * s012 * theta));

        c5[2] = DerivativeTerm(
            0., 2 * (qperpw * qperpy * s000 - qperpx * qperpz * s000 +
                q0y * q0z * s010 - qperpw * qperpx * s010 -
                qperpy * qperpz * s010 - q0y * q0y * s020 +
                qperpx * qperpx * s020 + qperpy * qperpy * s020 +
                q0x * q0z * (s000 - s100) - qperpw * qperpy * s100 +
                qperpx * qperpz * s100 +
                q0w * (q0y * (-s000 + s100) + q0x * (s010 - s110)) -
                q0y * q0z * s110 + qperpw * qperpx * s110 +
                qperpy * qperpz * s110 + q0y * q0y * s120 -
                qperpx * qperpx * s120 - qperpy * qperpy * s120 +
                q0x * q0x * (-s020 + s120)) *
            theta,
            2 * (qperpw * qperpy * s001 - qperpx * qperpz * s001 +
                q0y * q0z * s011 - qperpw * qperpx * s011 -
                qperpy * qperpz * s011 - q0y * q0y * s021 +
                qperpx * qperpx * s021 + qperpy * qperpy * s021 +
                q0x * q0z * (s001 - s101) - qperpw * qperpy * s101 +
                qperpx * qperpz * s101 +
                q0w * (q0y * (-s001 + s101) + q0x * (s011 - s111)) -
                q0y * q0z * s111 + qperpw * qperpx * s111 +
                qperpy * qperpz * s111 + q0y * q0y * s121 -
                qperpx * qperpx * s121 - qperpy * qperpy * s121 +
                q0x * q0x * (-s021 + s121)) *
            theta,
            2 * (qperpw * qperpy * s002 - qperpx * qperpz * s002 +
                q0y * q0z * s012 - qperpw * qperpx * s012 -
                qperpy * qperpz * s012 - q0y * q0y * s022 +
                qperpx * qperpx * s022 + qperpy * qperpy * s022 +
                q0x * q0z * (s002 - s102) - qperpw * qperpy * s102 +
                qperpx * qperpz * s102 +
                q0w * (q0y * (-s002 + s102) + q0x * (s012 - s112)) -
                q0y * q0z * s112 + qperpw * qperpx * s112 +
                qperpy * qperpz * s112 + q0y * q0y * s122 -
                qperpx * qperpx * s122 - qperpy * qperpy * s122 +
                q0x * q0x * (-s022 + s122)) *
            theta);
    }
}
void AnimatedTransform::Decompose(const Matrix4x4& m, Vector3f* T,
	Quaternion* Rquat, Matrix4x4* S)
{
	// Extract translation T from transformation matrix 
	T->x() = m.row(0)[3];
	T->y() = m.row(1)[3];
	T->z() = m.row(2)[3];

	// Compute new transformation matrix M without translation
	Matrix4x4 M = m;
	for (int i = 0; i < 3; ++i)
		M.row(i)[3] = M.row(3)[i] = 0.0;

	M.row(3)[3] = 1.0;

	// Extract rotation R from transformation matrix
	Float norm;
	int count = 0;
	Matrix4x4 R = M;
	do {
		// Compute next matrix Rnext in series 
		Matrix4x4 Rnext;
		Matrix4x4 Rit = (R).transpose().inverse();
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				Rnext.row(i)[j] = 0.5f * (R.row(i)[j] + Rit.row(i)[j]);

		// Compute norm of difference between R and Rnext 105
		norm = 0;
		for (int i = 0; i < 3; ++i) {
			Float n = std::abs(R.row(i)[0] - Rnext.row(i)[0]) +
				std::abs(R.row(i)[1] - Rnext.row(i)[1]) +
				std::abs(R.row(i)[2] - Rnext.row(i)[2]);
			norm = std::max(norm, n);
		}

		R = Rnext;
	} while (++count < 100 && norm > .0001);
	*Rquat = Quaternion(R);

	// Compute scale S using rotation and original matrix
	*S = R.inverse() * M;

	// Flip R[1] if needed to select shortest path
	if (R.row(0).dot(R.row(1)) < 0)
		R.row(0) = -R.row(0);
}

void AnimatedTransform::Interpolate(Float time, Transform* t) const
{
	// Handle boundary conditions for matrix interpolation 106
	if (!actuallyAnimated || time <= startTime) 
	{
		*t = *startTransform;
		return;
	}
	if (time >= endTime) 
	{
		*t = *endTransform;
		return;
	}

	Float dt = (time - startTime) / (endTime - startTime);

	// Interpolate translation at dt 106
	Vector3f trans = (1 - dt) * T[0] + dt * T[1];

	// Interpolate rotation at dt 107
	Quaternion rotate = Slerp(dt, R[0], R[1]);

	// Interpolate scale at dt 107
	Matrix4x4 scale;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			scale.row(i)[j] = Lerp(dt, S[0].row(i)[j], S[1].row(i)[j]);

	// Compute interpolated matrix as product of interpolated components 107
	*t = Translate(trans) * rotate.ToTransform() * Transform(scale);
}


Ray AnimatedTransform::operator()(const Ray& r) const 
{
    if (!actuallyAnimated || r.time <= startTime)
        return (*startTransform)(r);
    else if (r.time >= endTime)
        return (*endTransform)(r);
    else
    {
        Transform t;
        Interpolate(r.time, &t);
        return t(r);
    }
}

RayDifferential AnimatedTransform::operator()(const RayDifferential& r) const
{
    if (!actuallyAnimated || r.time <= startTime)
        return (*startTransform)(r);
    else if (r.time >= endTime)
        return (*endTransform)(r);
    else 
    {
        Transform t;
        Interpolate(r.time, &t);
        return t(r);
    }
}

Point3f AnimatedTransform::operator()(Float time, const Point3f& p) const
{
    if (!actuallyAnimated || time <= startTime)
        return (*startTransform)(p);
    else if (time >= endTime)
        return (*endTransform)(p);
    Transform t;
    Interpolate(time, &t);
    return t(p);
}

Vector3f AnimatedTransform::operator()(Float time, const Vector3f& v) const
{
    if (!actuallyAnimated || time <= startTime)
        return (*startTransform)(v);
    else if (time >= endTime)
        return (*endTransform)(v);
    Transform t;
    Interpolate(time, &t);
    return t(v);
}


Bounds3f AnimatedTransform::MotionBounds(const Bounds3f& b) const
{
	if (!actuallyAnimated)
		return (*startTransform)(b);
	if (hasRotation == false)
		return Union((*startTransform)(b), (*endTransform)(b));

	// Return motion bounds accounting for animated rotation 108
	Bounds3f bounds;
	for (int corner = 0; corner < 8; ++corner)
		bounds = Union(bounds, BoundPointMotion(b.Corner(corner)));
	return bounds;

}

Bounds3f AnimatedTransform::BoundPointMotion(const Point3f& p) const
{
	Bounds3f bounds((*startTransform)(p), (*endTransform)(p));
	Float cosTheta = Dot(R[0], R[1]);
	Float theta = std::acos(Clamp(cosTheta, -1, 1));
	for (int c = 0; c < 3; ++c)
	{
		Float zeros[4];
		int nZeros = 0;
		IntervalFindZeros(c1[c].Eval(p), c2[c].Eval(p), c3[c].Eval(p),
			c4[c].Eval(p), c5[c].Eval(p), theta,
			Interval(0., 1.), zeros, &nZeros);

		for (int i = 0; i < nZeros; ++i) 
		{
			Point3f pz = (*this)(Lerp(zeros[i], startTime, endTime), p);
			bounds = Union(bounds, pz);
		}
	}
	return bounds;
}
