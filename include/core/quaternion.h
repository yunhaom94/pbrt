#pragma once
#include "pbrt.h"

class Quaternion
{
public:
	Vector3f v;
	Float w;

public:
	Quaternion() : v(0, 0, 0), w(1) { }
    Quaternion(const Transform& t);
	~Quaternion() {}

    Quaternion& operator+=(const Quaternion& q) 
    {
        v += q.v;
        w += q.w;
        return *this;
    }

    friend Quaternion operator+(const Quaternion& q1, const Quaternion& q2)
    {
        Quaternion ret = q1;
        return ret += q2;
    }

    Quaternion& operator-=(const Quaternion& q) 
    {
        v -= q.v;
        w -= q.w;
        return *this;
    }
    Quaternion operator-() const
    {
        Quaternion ret;
        ret.v = -v;
        ret.w = -w;
        return ret;
    }

    friend Quaternion operator-(const Quaternion& q1, const Quaternion& q2) 
    {
        Quaternion ret = q1;
        return ret -= q2;
    }

    Quaternion& operator*=(Float f) 
    {
        v *= f;
        w *= f;
        return *this;
    }

    Quaternion operator*(Float f) const
    {
        Quaternion ret = *this;
        ret.v *= f;
        ret.w *= f;
        return ret;
    }

    Quaternion& operator/=(Float f) 
    {
        v /= f;
        w /= f;
        return *this;
    }

    Quaternion operator/(Float f) const
    {
        Quaternion ret = *this;
        ret.v /= f;
        ret.w /= f;
        return ret;
    }

    


    Transform ToTransform() const;
};

inline Quaternion operator*(Float f, const Quaternion& q) { return q * f; }


inline Float Dot(const Quaternion& q1, const Quaternion& q2) 
{
    return q1.v.dot(q2.v) + q1.w * q2.w;
}

inline Quaternion Normalize(const Quaternion & q) 
{
    return q / std::sqrt(Dot(q, q));
}


Quaternion Slerp(Float t, const Quaternion& q1,
    const Quaternion& q2);