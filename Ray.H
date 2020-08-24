#ifndef RAY_H
#define RAY_H
#include "slVector.H"
class Ray
{
public:
    SlVector3 e;
    SlVector3 d;
    // SlVector3 inv_d;
    int depth;
    Ray(const SlVector3 &_e, const SlVector3 &_d, int _depth = 0) : e(_e), d(_d), depth(_depth){};
};
#endif