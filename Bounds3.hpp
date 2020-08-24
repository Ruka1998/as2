#ifndef RAYTRACING_BOUNDS3_H
#define RAYTRACING_BOUNDS3_H
#include "Ray.H"
#include <limits>
#include <array>

class Bounds3
{
public:
    SlVector3 pMin, pMax; // two points to specify the bounding box
    Bounds3()
    {
        double minNum = std::numeric_limits<double>::lowest();
        double maxNum = std::numeric_limits<double>::max();
        pMax = SlVector3(minNum, minNum, minNum);
        pMin = SlVector3(maxNum, maxNum, maxNum);
    }
    Bounds3(const SlVector3 p) : pMin(p), pMax(p) {}
    Bounds3(const SlVector3 p1, const SlVector3 p2)
    {
        pMin = SlVector3(fmin(p1[0], p2[0]), fmin(p1[1], p2[1]), fmin(p1[2], p2[2]));
        pMax = SlVector3(fmax(p1[0], p2[0]), fmax(p1[1], p2[1]), fmax(p1[2], p2[2]));
    }

    SlVector3 Diagonal() const { return pMax - pMin; }
    int maxExtent() const
    {
        SlVector3 d = Diagonal();
        if (d[0] > d[1] && d[0] > d[2])
            return 0;
        else if (d[1] > d[2])
            return 1;
        else
            return 2;
    }

    double SurfaceArea() const
    {
        SlVector3 d = Diagonal();
        return 2 * (d[0] * d[1] + d[0] * d[2] + d[1] * d[2]);
    }

    SlVector3 Centroid() { return 0.5 * pMin + 0.5 * pMax; }
    Bounds3 Intersect(const Bounds3 &b)
    {
        return Bounds3(SlVector3(fmax(pMin[0], b.pMin[0]), fmax(pMin[1], b.pMin[1]),
                                 fmax(pMin[2], b.pMin[2])),
                       SlVector3(fmin(pMax[0], b.pMax[0]), fmin(pMax[1], b.pMax[1]),
                                 fmin(pMax[2], b.pMax[2])));
    }

    SlVector3 Offset(const SlVector3 &p) const
    {
        SlVector3 o = p - pMin;
        if (pMax[0] > pMin[0])
            o[0] /= pMax[0] - pMin[0];
        if (pMax[1] > pMin[1])
            o[1] /= pMax[1] - pMin[1];
        if (pMax[2] > pMin[2])
            o[2] /= pMax[2] - pMin[2];
        return o;
    }

    bool Overlaps(const Bounds3 &b1, const Bounds3 &b2)
    {
        bool x = (b1.pMax[0] >= b2.pMin[0]) && (b1.pMin[0] <= b2.pMax[0]);
        bool y = (b1.pMax[1] >= b2.pMin[1]) && (b1.pMin[1] <= b2.pMax[1]);
        bool z = (b1.pMax[2] >= b2.pMin[2]) && (b1.pMin[2] <= b2.pMax[2]);
        return (x && y && z);
    }

    bool Inside(const SlVector3 &p, const Bounds3 &b)
    {
        return (p[0] >= b.pMin[0] && p[0] <= b.pMax[0] && p[1] >= b.pMin[1] &&
                p[1] <= b.pMax[1] && p[2] >= b.pMin[2] && p[2] <= b.pMax[2]);
    }
    inline const SlVector3 &operator[](int i) const
    {
        return (i == 0) ? pMin : pMax;
    }

    inline bool IntersectP(const Ray &ray) const;
};

inline bool Bounds3::IntersectP(const Ray &ray) const
// inline bool Bounds3::IntersectP(const Ray &ray, const SlVector3 &invDir,
//                                 const std::array<int, 3> &dirisNeg) const
{
    // invDir: ray direction(x,y,z), invDir=(1.0/x,1.0/y,1.0/z), use this because Multiply is faster that Division
    // dirIsNeg: ray direction(x,y,z), dirIsNeg=[int(x>0),int(y>0),int(z>0)], use this to simplify your logic
    // TODO test if ray bound intersects
    SlVector3 invDir = 1.0 / ray.d;
    float tmin = (pMin[0] - ray.e[0]) * invDir[0];
    float tmax = (pMax[0] - ray.e[0]) * invDir[0];

    if (tmin > tmax)
    {
        std::swap(tmin, tmax);
    }

    float tymin = (pMin[1] - ray.e[1]) * invDir[1];
    float tymax = (pMax[1] - ray.e[1]) * invDir[1];

    if (tymin > tymax)
    {
        std::swap(tymin, tymax);
    }

    if ((tmin > tymax) || (tymin > tmax))
    {
        return false;
    }

    tmin = std::fmax(tymin, tmin);
    tmax = std::fmin(tymax, tmax);

    float tzmin = (pMin[2] - ray.e[2]) * invDir[2];
    float tzmax = (pMax[2] - ray.e[2]) * invDir[2];

    if (tzmin > tzmax)
    {
        std::swap(tzmin, tzmax);
    }

    if ((tmin > tzmax) || (tzmin > tmax))
    {
        return false;
    }

    tmin = std::fmax(tzmin, tmin);
    tmax = std::fmin(tzmax, tmax);

    return true;
}

inline Bounds3 Union(const Bounds3 &b1, const Bounds3 &b2)
{
    Bounds3 ret;
    ret.pMin = min(b1.pMin, b2.pMin);
    ret.pMax = max(b1.pMax, b2.pMax);
    return ret;
}

inline Bounds3 Union(const Bounds3 &b, const SlVector3 &p)
{
    Bounds3 ret;
    ret.pMin = min(b.pMin, p);
    ret.pMax = max(b.pMax, p);
    return ret;
}

#endif // RAYTRACING_BOUNDS3_H
