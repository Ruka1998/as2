#include "Object.H"

#define TVALID(t) (t >= t0 && t < t1)

bool Triangle::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const
{
    // Step 1 Ray-triangle test
    SlVector3 E1 = b - a, E2 = c - a;
    SlVector3 S = r.e - a, S1 = cross(r.d, E2), S2 = cross(S, E1);
    double denominator = dot(S1, E1);
    if (denominator < 1e-10 && denominator > -1e-10)
        return false;
    hr.t = dot(S2, E2) / denominator;
    if (!TVALID(hr.t))
        return false;
    hr.beta = dot(S1, S) / denominator;
    if (hr.beta < 0 || hr.beta > 1)
        return false;
    hr.gamma = dot(S2, r.d) / denominator;
    if (hr.gamma < 0 || hr.gamma > 1)
        return false;
    hr.alpha = 1 - hr.beta - hr.gamma;
    if (hr.alpha < 0)
        return false;
    hr.p = r.e + hr.t * r.d;
    hr.v = r.e - hr.p;
    hr.n = n;
    normalize(hr.v);
    return true;
}

bool TrianglePatch::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const
{
    bool temp = Triangle::intersect(r, t0, t1, hr);
    if (temp)
    {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
    }
    return temp;
}

bool Sphere::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const
{
    // Step 1 Sphere-triangle test
    SlVector3 tmp1 = r.e - c;
    double tmp2 = dot(r.d, tmp1);
    double tmp3 = dot(tmp1, tmp1) - rad * rad;
    double delta = tmp2 * tmp2 - tmp3;

    if (delta < 0)
        return false;
    double tmp4 = dot(r.d, r.d);
    if (delta == 0)
    {
        hr.t = -tmp2 / tmp4;
        if (hr.t < t0 || hr.t >= t1)
            return false;
        hr.p = r.e + r.d * hr.t;
        hr.n = hr.p - c;
        normalize(hr.n);
        hr.v = r.e - hr.p;
        normalize(hr.v);
        return true;
    }
    else
    {
        double tmp5 = sqrt(delta);
        double tx1 = (-tmp2 - tmp5) / tmp4, tx2 = (-tmp2 + tmp5) / tmp4;
        if (!TVALID(tx1) && !TVALID(tx2))
            return false;
        else if (TVALID(tx1))
            hr.t = tx1;
        else if (TVALID(tx2))
            hr.t = tx2;
        else
            hr.t = std::min(tx1, tx2);
        hr.p = r.e + r.d * hr.t;
        hr.n = hr.p - c;
        normalize(hr.n);
        hr.v = r.e - hr.p;
        normalize(hr.v);
        return true;
    }
}