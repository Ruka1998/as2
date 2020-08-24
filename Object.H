#ifndef OBJECT_H
#define OBJECT_H
#include "Bounds3.hpp"

class Fill
{
public:
    SlVector3 color;
    double kd, ks, shine, t, ior;
};

class Light
{
public:
    SlVector3 p, c;
};

class HitRecord
{
public:
    double t, alpha, beta, gamma;
    SlVector3 p, n, v;
    Fill f;
    int raydepth;
};

class Surface
{
public:
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const = 0;
    virtual ~Surface(){};
    virtual Bounds3 getBounds() = 0;
};

class Triangle : public Surface
{
    SlVector3 a, b, c;
    SlVector3 n;

public:
    Triangle(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c) : a(_a), b(_b), c(_c){
        SlVector3 e1 = b-a, e2 = c-a;
        normalize(e1);
        normalize(e2);
        n = cross(e1,e2);
    };
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    Bounds3 getBounds()
    {
        return Union(Bounds3(a, b), c);
    }
};

class TrianglePatch : public Triangle
{
    SlVector3 n1, n2, n3;

public:
    TrianglePatch(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c,
                  const SlVector3 &_n1, const SlVector3 &_n2, const SlVector3 &_n3)
        : Triangle(_a, _b, _c), n1(_n1), n2(_n2), n3(_n3){};
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
};

class Sphere : public Surface
{
    SlVector3 c;
    double rad;

public:
    Sphere(const SlVector3 &_c, double _r) : c(_c), rad(_r){};
    bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    Bounds3 getBounds()
    {
        return Bounds3(SlVector3(c[0] - rad, c[1] - rad, c[2] - rad),
                       SlVector3(c[0] + rad, c[1] + rad, c[2] + rad));
    }
};
#endif