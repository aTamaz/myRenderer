// smallpt, a Path Tracer by Kevin Beason, 2008
// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
//        Remove "-fopenmp" for g++ version < 4.2
// Usage: time ./smallpt 5000 && xv image.ppm
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include "lodepng.h"

using namespace std;

void PressEnterToContinue()
{
    int c;
    printf("Press ENTER to continue... ");
    fflush(stdout);
    do c = getchar(); while ((c != '\n') && (c != EOF));
}

//double erand48(const unsigned short &Xi) {
    //srand(Xi);
//    return rand() / (double)RAND_MAX;
//}

struct Vec {
    double x, y, z;                  // position, also color (r,g,b)
    //Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
    Vec operator+(const Vec &b) const {
        return Vec(x + b.x, y + b.y, z + b.z);
    }
    Vec operator-(const Vec &b) const {
        return Vec(x - b.x, y - b.y, z - b.z);
    }
    Vec operator*(double b) const {
        return Vec(x*b, y*b, z*b);
    }
    Vec mult(const Vec &b) const {
        return Vec(x * b.x, y * b.y, z * b.z);
    }
    Vec & norm() {
        return *this = *this * (1/sqrt(x*x + y*y + z*z));
    }
    double dot(const Vec &b) const { // cross product
        return x*b.x + y*b.y + z*b.z;
    }
    Vec operator%(Vec &b) {
        return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
    }
};
struct Ray {
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
    double rad;       // radius
    Vec p, e, c;      // position, emission, color
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
        rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double eps = 1e-4;
        double b = op.dot(r.d);
        double det = b * b - op.dot(op) + rad * rad;
        if (det < 0) {
            return 0;
        } else {
            det = sqrt(det);
        }
        double t = b - det;
        if (t > eps) return t;
        t = b + det;
        if (t > eps) return t;
        return 0;
    }
};
// x horizontal: left2right
// y vertical: bottom2up
// z depth: back2front (out)
Sphere spheres[] = { //Scene: radius, position, emission, color, material
    Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF), // Left
    Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF), // Right
    Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF), // Back
    Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF), // Front
    Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF), // Bottom
    Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF), // Top
    Sphere(16.5,Vec(27,16.5,90),       Vec(),Vec(1,1,1)*.999, SPEC), // Mirror
    Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR), // Glas
//    Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF)  // Light
    Sphere(1.5, Vec(50,81.6-16.5,81.6),Vec(4,4,4)*100,  Vec(), DIFF),  // small light
};
int numSpheres = sizeof(spheres)/sizeof(Sphere);
inline double clamp(double x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}
inline int toInt(double x)
{
    return int(pow(clamp(x),1/2.2)*255+.5);
}
void saveImage(Vec* image, int width, int height, int numSamples, int seconds)
{
    vector<unsigned char> pngImage;
    pngImage.resize(width*height*4);
    for (int i = 0; i < width*height; i++) {
        pngImage[4*i + 0] = toInt(image[i].x);
        pngImage[4*i + 1] = toInt(image[i].y);
        pngImage[4*i + 2] = toInt(image[i].z);
        pngImage[4*i + 3] = 255;
    }

    vector<unsigned char> png;
    unsigned error = lodepng::encode(png, pngImage, width, height);
    if (error) {
        cout << "encoder error: " << error << ": " << lodepng_error_text(error) << endl;
        return;
    }

    time_t t = time(0);
    struct tm *now = localtime(&t);
    ostringstream ost;
    //ost << "/Users/tamaz/Desktop/"
    ost << "../results/"
        << (now->tm_year + 1900)                            << "."
        << setw(2) << setfill('0') << now->tm_mon + 1       << "."
        << setw(2) << setfill('0') << now->tm_mday          << "_"
        << setw(2) << setfill('0') << now->tm_hour          << "."
        << setw(2) << setfill('0') << now->tm_min           << "."
        << setw(2) << setfill('0') << now->tm_sec           << "_"
        << numSamples << "-spp_" << seconds << "-secs.png";

    lodepng::save_file(png, ost.str());
}
inline bool intersect(const Ray &r, double &t, int &id)
{
    double inf = t = 1e20;
    for (int i = 0; i < numSpheres; i++) {
        double d = spheres[i].intersect(r);
        if (d != 0 && d < t) {
            t = d;
            id = i;
        }
    }
    return t<inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi, int E = 1)
{
    double t;                               // distance to intersection
    int id(0);                              // id of intersected object
    if (!intersect(r, t, id)) return Vec(); // if miss, return black
    const Sphere &obj = spheres[id];        // the hit object
    Vec x = r.o + r.d * t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    Vec f = obj.c;
    double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
    if (++depth > 5 || !p) {
        double rnd = erand48(Xi); // Russian roulette
        if (rnd < p && depth < 100) {
            f = f*(1/p);
        } else {
            return obj.e*E;
        }
    }
    if (obj.refl == DIFF) {                  // Ideal DIFFUSE reflection
        double r1 = 2*M_PI*erand48(Xi);
        double r2 = erand48(Xi);
        double r2s = sqrt(r2);
        Vec w = nl;
        Vec u = ((fabs(w.x) > .1 ? Vec(0,1) : Vec(1)) % w).norm();
        Vec v = w % u;
        Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
        // Loop over any lights
        Vec e;
        for (int i=0; i<numSpheres; i++){
            const Sphere &s = spheres[i];
            if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0) continue; // skip non-lights

            Vec sw = s.p - x;
            Vec su = ((fabs(sw.x) > .1 ? Vec(0,1) : Vec(1)) % sw).norm();
            Vec sv = sw % su;
            double cos_a_max = sqrt(1 - s.rad*s.rad/(x - s.p).dot(x - s.p));
            double eps1 = erand48(Xi);
            double eps2 = erand48(Xi);
            double cos_a = 1 - eps1 + eps1 * cos_a_max;
            double sin_a = sqrt(1 - cos_a * cos_a);
            double phi = 2*M_PI*eps2;
            Vec l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
            l.norm();
            if (intersect(Ray(x,l), t, id) && id == i){  // shadow ray
                double omega = 2*M_PI*(1 - cos_a_max);
                e = e + f.mult(s.e*l.dot(nl)*omega)*M_1_PI;  // 1/pi for brdf
            }
        }
        return obj.e*E + f.mult(radiance(Ray(x,d),depth,Xi, 0));
    }
    if (obj.refl == SPEC)                       // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
    Ray reflRay(x, r.d - n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1;
    double nt = 1.5;
    double nnt = into ? nc/nt : nt/nc;
    double ddn = r.d.dot(nl);
    double cos2t = 1 - nnt*nnt*(1-ddn*ddn);
    if (cos2t < 0)    // Total internal reflection
        return obj.e + f.mult(radiance(reflRay,depth,Xi));
    Vec tdir = (r.d*nnt - n * ((into?1:-1) * (ddn*nnt+sqrt(cos2t)))).norm();
    double a = nt-nc, b = nt+nc, R0 = a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
    double Re = R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
    return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
        radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
        radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
void render(int numSamples, int w, int h)
{
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
    Vec cx = Vec(w * .5135/h), cy = (cx % cam.d).norm() * .5135;
    Vec r;
    Vec *c = new Vec[w * h];
    //clock_t begin = clock();
    int t0 = time(NULL);
    //omp_set_dynamic(0);
    #pragma omp parallel for schedule(dynamic, 1) private(r) num_threads(4) // OpenMP
    for (int y = 0; y < h; y++) { // Loop over image rows
        if (y % 10 == 0)
            fprintf(stderr, "\rRendering (%d spp) %5.2f%%", numSamples*4, 100.*y/(h-1));
        unsigned short Xi[3] = {0, 0, y*y*y};
        //srand(Xi);
        for (int x = 0; x < w; x++) {  // Loop cols
            int i = (h - y - 1) * w + x;
            for (int sy = 0; sy < 2; sy++) {     // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++) {        // 2x2 subpixel cols
                    r = Vec();
                    for (int s = 0; s < numSamples; s++) {
                        double r1 = 2*erand48(Xi);
                        double dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                        double r2 = 2*erand48(Xi);
                        double dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
                        Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                                cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./numSamples);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
            }
        }
    }
    //clock_t end = clock();
    //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    int elapsed_secs = time(NULL) - t0;
    saveImage(c, w, h, numSamples * 4, elapsed_secs);
    delete[] c;
    cout << endl << "Done in " << elapsed_secs << " seconds!" << endl;
}
int main(int argc, char *argv[])
{
    int defaultNumSamples = 400;
    int numSamples = argc==2 ? atoi(argv[1]) : defaultNumSamples;  // # samples per subpixel
    numSamples /= 4;
    int w = 1024, h = 768;
    render(numSamples, w, h);
    return 0;
}
