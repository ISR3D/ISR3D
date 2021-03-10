#ifndef KDTREE_HIT_H
#define KDTREE_HIT_H

#include <cstddef>
#include <iostream>
#include <functional>

namespace kdtree {

/// Data structure used to store the index of an object matching a spatial query, along with
/// its distance squared to the query object.
template<typename float_t>
struct Hit {
    size_t dIdx;        // direct index into data
    float_t distSqr;    // distance squared

    Hit(size_t dIdx_, float_t distSqr_) : dIdx(dIdx_), distSqr(distSqr_) { }
    // note: we do not guard against calling the constructor with swapped arguments,
    // which is possible due to implicit type conversion; be careful.
};

template<typename float_t>
bool operator<(Hit<float_t> const& a, Hit<float_t> const& b) {
    return ((a.distSqr) < (b.distSqr));
}

template<typename float_t>
bool operator==(Hit<float_t> const& a, Hit<float_t> const& b) {
    return (a.dIdx==b.dIdx && a.distSqr==b.distSqr);
}

template<typename float_t>
std::ostream & operator<<(std::ostream & ostr, Hit<float_t> result)
{
    ostr << result.dIdx << ":" << result.distSqr << "  ";
    return ostr;
}

} /* namespace kdtree */

#endif
