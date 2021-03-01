//
//  Simple 3D vector class
//
#pragma once

struct Vector {
    double  x, y, z;

            Vector() {};
            Vector(double a, double b, double c) : x(a), y(b), z(c) {};
    Vector operator +(Vector v) { return Vector(x+v.x, y+v.y, z+v.z); }
    Vector operator -(Vector v) { return Vector(x-v.x, y-v.y, z-v.z); }
    Vector operator *(double f) { return Vector(x*f, y*f, z*f); }
    double abs() { return sqrt(x*x + y*y + z*z); }
};

inline double
DotProduct(const Vector &vec1, const Vector &vec2)
{
    return  vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
}

inline Vector
CrossProduct(const Vector &vec1, const Vector &vec2)
{
    return Vector(vec1.y*vec2.z-vec1.z*vec2.y,
                  vec1.z*vec2.x-vec1.x*vec2.z,
                  vec1.x*vec2.y-vec1.y*vec2.x);
}
