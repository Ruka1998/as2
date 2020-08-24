#include "trace.H"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
#include <values.h>
#include <chrono>
#define MAX DBL_MAX
#endif

#define MAX_DEPTH 4

// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3 &a, const SlVector3 &b, const SlVector3 &c)
{
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
           b[0] * (c[1] * a[2] - a[1] * c[2]) +
           c[0] * (a[1] * b[2] - b[1] * a[2]);
}

inline double sqr(double x) { return x * x; }

SlVector3 reflect(SlVector3 i, SlVector3 n)
{
    normalize(i);
    normalize(n);
    return 2 * dot(i, n) * n - i;
}

SlVector3 refract(SlVector3 i, SlVector3 n, const float ior)
{
    normalize(i);
    normalize(n);
    float cosi = dot(i, n);
    float etai = 1, etat = ior;
    if (cosi < 0)
    {
        cosi = -cosi;
    }
    else
    {
        std::swap(etai, etat);
        n = -n;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? SlVector3(0.0) : eta * i + (eta * cosi - sqrtf(k)) * n;
}

float fresnel(SlVector3 I, SlVector3 N, const float &ior)
{
    float kr;
    normalize(I);
    normalize(N);
    float cosi = dot(I, N);
    float etai = 1, etat = ior;
    if (cosi > 0)
    {
        std::swap(etai, etat);
    }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1)
    {
        kr = 1;
    }
    else
    {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
    return kr;
}

Tracer::Tracer(const std::string &fname)
{
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in)
    {
        getline(in, line);
        switch (line[0])
        {
        case 'b':
        {
            std::stringstream ss(line);
            ss >> ch >> bcolor[0] >> bcolor[1] >> bcolor[2];
            break;
        }

        case 'v':
        {
            getline(in, line);
            std::string junk;
            std::stringstream fromss(line);
            fromss >> junk >> eye[0] >> eye[1] >> eye[2];

            getline(in, line);
            std::stringstream atss(line);
            atss >> junk >> at[0] >> at[1] >> at[2];

            getline(in, line);
            std::stringstream upss(line);
            upss >> junk >> up[0] >> up[1] >> up[2];

            getline(in, line);
            std::stringstream angless(line);
            angless >> junk >> angle;

            getline(in, line);
            std::stringstream hitherss(line);
            hitherss >> junk >> hither;

            getline(in, line);
            std::stringstream resolutionss(line);
            resolutionss >> junk >> res[0] >> res[1];
            break;
        }

        case 'p':
        {
            bool patch = false;
            std::stringstream ssn(line);
            unsigned int nverts;
            if (line[1] == 'p')
            {
                patch = true;
                ssn >> ch;
            }
            ssn >> ch >> nverts;
            std::vector<SlVector3> vertices;
            std::vector<SlVector3> normals;
            for (unsigned int i = 0; i < nverts; i++)
            {
                getline(in, line);
                std::stringstream ss(line);
                SlVector3 v, n;
                if (patch)
                    ss >> v[0] >> v[1] >> v[2] >> n[0] >> n[1] >> n[2];
                else
                    ss >> v[0] >> v[1] >> v[2];
                vertices.push_back(v);
                normals.push_back(n);
            }
            bool makeTriangles = false;
            if (vertices.size() == 3)
            {
                if (patch)
                {
                    surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                                                                                    normals[0], normals[1], normals[2]),
                                                                  fill));
                }
                else
                {
                    surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                }
            }
            else if (vertices.size() == 4)
            {
                SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
                SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
                SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
                SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
                if (dot(n0, n1) > 0 && dot(n0, n2) > 0 && dot(n0, n3) > 0)
                {
                    makeTriangles = true;
                    if (patch)
                    {
                        surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                                                                                        normals[0], normals[1], normals[2]),
                                                                      fill));
                        surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3],
                                                                                        normals[0], normals[2], normals[3]),
                                                                      fill));
                    }
                    else
                    {
                        surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                        surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
                    }
                }
                if (!makeTriangles)
                {
                    std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
                }
            }
            break;
        }

        case 's':
        {
            std::stringstream ss(line);
            SlVector3 c;
            double r;
            ss >> ch >> c[0] >> c[1] >> c[2] >> r;
            surfaces.push_back(std::pair<Surface *, Fill>(new Sphere(c, r), fill));
            break;
        }

        case 'f':
        {
            std::stringstream ss(line);
            ss >> ch >> fill.color[0] >> fill.color[1] >> fill.color[2] >> fill.kd >> fill.ks >> fill.shine >> fill.t >> fill.ior;
            break;
        }

        case 'l':
        {
            std::stringstream ss(line);
            Light l;
            ss >> ch >> l.p[0] >> l.p[1] >> l.p[2];
            if (!ss.eof())
            {
                ss >> l.c[0] >> l.c[1] >> l.c[2];
                coloredlights = true;
            }
            lights.push_back(l);
            break;
        }

        default:
            break;
        }
    }
    if (!coloredlights)
        for (unsigned int i = 0; i < lights.size(); i++)
            lights[i].c = 1.0 / sqrt(lights.size());
    im = new SlVector3[res[0] * res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
    for (int i = 0; i < surfaces.size(); i++)
    {
        objects.push_back(&surfaces[i]);
    }
}

Tracer::~Tracer()
{
    if (im)
        delete[] im;
    for (unsigned int i = 0; i < surfaces.size(); i++)
        delete surfaces[i].first;
}

SlVector3 Tracer::shade(const HitRecord &hr, double t0, double t1) const
{
    if (color)
        return hr.f.color;

    SlVector3 color(0.0);
    HitRecord dummy;

    for (unsigned int i = 0; i < lights.size(); i++)
    {
        const Light &light = lights[i];
        bool shadow = false;

        // Step 3 Check for shadows here
        HitRecord tmpHR;
        SlVector3 PtoL = light.p - hr.p;
        normalize(PtoL);
        Ray shadowRay(hr.p, PtoL);
        if (bvh->Intersect(shadowRay, 0.005, MAX, tmpHR))
        {
            shadow = true;
        }

        if (!shadow)
        {
            // Step 2 do shading here
            SlVector3 L = light.p - hr.p;
            normalize(L);
            SlVector3 V = hr.v;
            SlVector3 H = L + V;
            normalize(H);
            SlVector3 diffuse = std::max(dot(L, hr.n), 0.0) * hr.f.kd;
            // SlVector3 specular = pow(dot(hr.n, H), hr.f.shine) * hr.f.ks;
            SlVector3 specular = pow(std::max(dot(hr.n, H), 0.0), hr.f.shine) * hr.f.ks;
            // color += specular;
            color += (diffuse * hr.f.color + specular) * light.c;
        }
    }

    // Step 4 Add code for computing reflection color here
    int kr = 1;

    Ray refl(hr.p, reflect(hr.v, hr.n), hr.raydepth + 1);
    SlVector3 reflectColor = trace(refl, t0, t1);
    // Step 5 Add code for computing refraction color here
    SlVector3 refractColor;
    if (hr.f.ior != 0.0)
    {
        Ray refr(hr.p, refract(-hr.v, hr.n, hr.f.ior), hr.raydepth + 1);
        refractColor = trace(refr, t0, t1);
        kr = fresnel(-hr.v, hr.n, hr.f.ior);
    }
    // reflectColor = SlVector3(0.0);
    return color + reflectColor * kr + refractColor * (1 - kr);
    // return color + reflectColor * (1 - hr.f.t) + refractColor * hr.f.t;
}

SlVector3 Tracer::trace(const Ray &r, double t0, double t1) const
{
    if (r.depth > MAX_DEPTH)
        return SlVector3(0.0, 0.0, 0.0);
    HitRecord hr;
    SlVector3 color(bcolor);

    bool hit = false;

    // Step 1 See what a ray hits
    hit = bvh->Intersect(r, t0, t1, hr);
    hr.raydepth = r.depth;

    if (hit)
        color = shade(hr, 0.00001, MAX);
    return color;
}

void Tracer::traceImage()
{
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up, w);
    normalize(u);
    SlVector3 v = cross(w, u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((M_PI / 180.0) * (angle / 2.0)) * d;
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3 *pixel = im;

    for (unsigned int j = 0; j < res[1]; j++)
    {
        for (unsigned int i = 0; i < res[0]; i++, pixel++)
        {

            SlVector3 result(0.0, 0.0, 0.0);

            for (int k = 0; k < samples; k++)
            {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r - l) * (i + rx) / res[0];
                double y = b + (t - b) * (j + ry) / res[1];
                SlVector3 dir = -d * w + x * u + y * v;

                Ray r(eye, dir);
                normalize(r.d);

                result += trace(r, hither, MAX);
            }
            (*pixel) = result / samples;
        }
        // std::cout << (int)(j*100 / res[1]) << "%" << std::endl;
    }
}

void Tracer::writeImage(const std::string &fname)
{
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out << "P6"
        << "\n"
        << res[0] << " " << res[1] << "\n"
        << 255 << "\n";
    SlVector3 *pixel = im;
    char val;
    for (unsigned int i = 0; i < res[0] * res[1]; i++, pixel++)
    {
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write(&val, sizeof(unsigned char));
    }
    out.close();
}

void Tracer::buildBVH()
{
    std::cout << "Generating BVH...\n";
    auto start = std::chrono::system_clock::now();
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds\n";
}

int main(int argc, char *argv[])
{
    int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    while ((c = getopt(argc, argv, "a:s:d:c")) != -1)
    {
        switch (c)
        {
        case 'a':
            aperture = atof(optarg);
            break;
        case 's':
            samples = atoi(optarg);
            break;
        case 'c':
            color = true;
            break;
        case 'd':
            maxraydepth = atoi(optarg);
            break;
        default:
            abort();
        }
    }

    // if (argc - optind != 2)
    // {
    //     std::cout << "usage: trace [opts] input.nff output.ppm" << std::endl;
    //     for (unsigned int i = 0; i < argc; i++)
    //         std::cout << argv[i] << std::endl;
    //     exit(0);
    // }
    Tracer tracer("InputFiles/refract.nff");

    // Tracer tracer(argv[optind++]);
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;
    // tracer.color = true;

    tracer.buildBVH();

    auto start = std::chrono::system_clock::now();
    tracer.traceImage();
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds\n";
    tracer.traceImage();
    tracer.writeImage("output2.ppm");
    // tracer.writeImage(argv[optind++]);
};