#include "core/quaternion.h"
#include "core/transform.h"

Transform Quaternion::ToTransform() const 
{
    Float xx = v.x() * v.x(), yy = v.y() * v.y(), zz = v.z() * v.z();
    Float xy = v.x() * v.y(), xz = v.x() * v.z(), yz = v.y() * v.z();
    Float wx = v.x() * w, wy = v.y() * w, wz = v.z() * w;

    Matrix4x4 m;
    m.row(0)[0] = 1 - 2 * (yy + zz);
    m.row(0)[1] = 2 * (xy + wz);
    m.row(0)[2] = 2 * (xz - wy);
    m.row(1)[0] = 2 * (xy - wz);
    m.row(1)[1] = 1 - 2 * (xx + zz);
    m.row(1)[2] = 2 * (yz + wx);
    m.row(2)[0] = 2 * (xz + wy);
    m.row(2)[1] = 2 * (yz - wx);
    m.row(2)[2] = 1 - 2 * (xx + yy);

    // Transpose since we are left-handed.  Ugh.
    return Transform(Transpose(m), m);
}

Quaternion::Quaternion(const Transform& t)
{
    const Matrix4x4& m = t.m;
    Float trace = m.row(0)[0] + m.row(1)[1] + m.row(2)[2];
    if (trace > 0.f)
    {
        // Compute w from matrix trace, then xyz
        // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
        Float s = std::sqrt(trace + 1.0f);
        w = s / 2.0f;
        s = 0.5f / s;
        v.x() = (m.row(2)[1] - m.row(1)[2]) * s;
        v.y() = (m.row(0)[2] - m.row(2)[0]) * s;
        v.z() = (m.row(1)[0] - m.row(0)[1]) * s;
    }
    else {
        // Compute largest of $x$, $y$, or $z$, then remaining components
        const int nxt[3] = { 1, 2, 0 };
        Float q[3];
        int i = 0;
        if (m.row(1)[1] > m.row(0)[0]) i = 1;
        if (m.row(2)[2] > m.row(i)[i]) i = 2;
        int j = nxt[i];
        int k = nxt[j];
        Float s = std::sqrt((m.row(i)[i] - (m.row(j)[j] + m.row(k)[k])) + 1.0f);
        q[i] = s * 0.5f;
        if (s != 0.f) s = 0.5f / s;
        w = (m.row(k)[j] - m.row(j)[k]) * s;
        q[j] = (m.row(j)[i] + m.row(i)[j]) * s;
        q[k] = (m.row(k)[i] + m.row(i)[k]) * s;
        v.x() = q[0];
        v.y() = q[1];
        v.z() = q[2];
    }
}

Quaternion Slerp(Float t, const Quaternion& q1,
    const Quaternion& q2) 
{
    Float cosTheta = Dot(q1, q2);
    if (cosTheta > 0.9995)
        return Normalize((1 - t) * q1 + t * q2);
    else 
    {
        Float theta = std::acos(Clamp(cosTheta, -1, 1));
        Float thetap = theta * t;
        Quaternion qperp = Normalize(q2 - q1 * cosTheta);
        return q1 * std::cos(thetap) + qperp * std::sin(thetap);
    }
}