#ifndef KDTREE_GEOMETRY_H
#define KDTREE_GEOMETRY_H

#include <cmath>
#include <cassert>
#include <iostream>

namespace kdtree {

/// A point in k-dimensional Cartesian space.
template<size_t k, typename float_t>
struct Point {
    float_t x[k];

    Point() {
        for (size_t i=0; i<k; i++) x[i]=float_t();
    }
    Point(Point const& rhs) {
        for (size_t i=0; i<k; i++) x[i]=rhs.x[i];
    }
    Point& operator=(Point const& rhs) {
        for (size_t i=0; i<k; i++) x[i]=rhs.x[i];
        return *this;
    }
    Point& operator+=(Point const& rhs) {
        for (size_t i=0; i<k; i++) x[i]+=rhs.x[i];
        return *this;
    }
    Point& operator-=(Point const& rhs) {
        for (size_t i=0; i<k; i++) x[i]-=rhs.x[i];
        return *this;
    }
    Point& operator*=(double const& rhs) {
        for (size_t i=0; i<k; i++) x[i]*=rhs;
        return *this;
    }
    Point& operator/=(double const& rhs) {
        for (size_t i=0; i<k; i++) x[i]/=rhs;
        return *this;
    }
    float_t& operator[](size_t i) {
        assert(i < k);
        return x[i];
    }
    float_t const& operator[](size_t i) const {
        assert(i < k);
        return x[i];
    }

    Point(float_t x0, float_t x1) {
        assert(k==2); // do not use this constructor for k!=2
        x[0]=x0;
        x[1]=x1;
    }
    Point(float_t x0, float_t x1, float_t x2) {
        assert(k==3); // do not use this constructor for k!=3
        x[0]=x0;
        x[1]=x1;
        x[2]=x2;
    }
};

template<size_t k, typename float_t>
bool operator==(Point<k, float_t> const& a, Point<k, float_t> const& b) {
    for (size_t i=0; i<k; i++) {
        if (a.x[i]!=b.x[i]) return false;
    }
    return true;
}

template<size_t k, typename float_t>
inline Point<k, float_t> operator-(Point<k, float_t> lhs, Point<k, float_t> const& rhs)
{
    lhs -= rhs;
    return lhs;
}

/// Dot product
template<size_t k, typename float_t>
inline float_t operator*(Point<k, float_t> const& lhs, Point<k, float_t> const& rhs)
{
    float_t temp = 0;
    for (size_t i=0; i<k; i++) {
        temp += lhs[i] * rhs [i];
    }
    return temp;
}

/// Cross product
template<typename float_t>
inline Point<3, float_t> cross(Point<3, float_t> lhs, Point<3, float_t> const& rhs)
{
    return Point<3, float_t> (lhs[1] * rhs [2] - lhs[2] * rhs [1],
                              lhs[2] * rhs [0] - lhs[0] * rhs [2],
                              lhs[0] * rhs [1] - lhs[1] * rhs [0]
                              );
}

template<size_t k, typename float_t>
inline Point<k, float_t> operator*(Point<k, float_t> const& lhs, float_t const& rhs)
{
    Point<k, float_t> temp = lhs;
    temp *= rhs;
    return temp;
}

template<size_t k, typename float_t>
inline Point<k, float_t> operator*(float_t const& lhs, Point<k, float_t> rhs)
{
    return rhs * lhs;
}

template<size_t k, typename float_t>
inline Point<k, float_t> operator/(Point<k, float_t> const& lhs, float_t const& rhs)
{
  Point<k, float_t> temp = lhs;
  temp /= rhs;
  return temp;
}

template<size_t k, typename float_t>
bool operator!=(Point<k, float_t> const& a, Point<k, float_t> const& b) {
    return !(a==b);
}

template<size_t k, typename float_t>
std::ostream& operator<<(std::ostream& ostr, Point<k, float_t> const& p)
{
    for (size_t i=0; i<k; i++) ostr << p.x[i] << " ";
    return ostr;
}

template<int nDim>
inline double normSqr(Point<nDim, double> const& v)
{
    double result = 0.0;
    for (int iDim=0; iDim<nDim; iDim++) {
        result += v[iDim]*v[iDim];
    }
    return result;
}

template<int nDim>
inline double norm(Point<nDim, double> const& v)
{
    return sqrt(normSqr<nDim>(v) );
}

/// A hyper-box in k-dimensional Cartesian space, represented by min and max points.
template<size_t k, typename float_t>
struct Box {
    typedef Point<k, float_t> point_t;

    point_t lo, hi;

    Box(point_t const& lo_, point_t const& hi_) : lo(lo_), hi(hi_) {
        for (size_t i=0; i<k; i++) { assert(lo.x[i] <= hi.x[i]); }
    }
};

template<size_t k, typename float_t>
bool operator==(Box<k, float_t> const& a, Box<k, float_t> const& b) {
    return (a.lo==b.lo && a.hi==b.hi);
}
template<size_t k, typename float_t>
bool operator!=(Box<k, float_t> const& a, Box<k, float_t> const& b) {
    return !(a==b);
}

template<size_t k, typename float_t>
std::ostream& operator<<(std::ostream& ostr, Box<k, float_t> const& box)
{
    ostr << box.lo << " " << box.hi << " ";
    return ostr;
}


} /* namespace kdtree */

#endif
