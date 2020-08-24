#ifndef TRACE_H
#define TRACE_H

#include "slVector.H"
#include "Ray.H"
#include "Object.H"
#include "Bounds3.hpp"
#include "BVH.hpp"
#include <vector>


class Tracer
{
    SlVector3 bcolor, eye, at, up;
    double angle, hither;
    unsigned int res[2];
    std::vector<std::pair<Surface *, Fill>> surfaces;
    std::vector<std::pair<Surface *, Fill>*> objects;
    std::vector<Light> lights;
    double shadowbias;

    SlVector3 *im;

public:
    Tracer(const std::string &fname);
    ~Tracer();
    void traceImage();
    SlVector3 trace(const Ray &ray, double t0, double t1) const;
    SlVector3 shade(const HitRecord &hr, double t0, double t1) const;
    void writeImage(const std::string &fname);

    bool color;
    int samples;
    double aperture;
    int maxraydepth;

    BVHAccel *bvh;
    void buildBVH();
};

#endif
