#ifndef KDTREE_EUCLIDEAN_METRICS_H
#define KDTREE_EUCLIDEAN_METRICS_H

#include <cmath>
#include <limits>
#include "geometry.h"

namespace kdtree {

/// Euclidean distance metrics for points and boxes in k-dimensional space.
template<size_t k, typename float_t>
class EuclideanMetrics {
public:
    typedef Point<k, float_t> point_t;
    typedef Box<k, float_t> box_t;

    /// Distance squared between two points.
    static float_t distSqr(point_t const& a, point_t const& b) {
        float_t dSqr = 0;
        for (size_t i=0; i<k; i++) dSqr += ((a.x[i]-b.x[i])*(a.x[i]-b.x[i]));
        return dSqr;
    }

    /// Distance between two points.
    static float_t dist(point_t const& a, point_t const& b) {
        return sqrt(distSqr(a,b) );
    }

    /// Shortest distance squared between a hyper-box and a point.
    static float_t distSqr(box_t const& box, point_t const& p) {
        float_t dSqr = 0;
        for (size_t i=0; i<k; i++) {
            if (p.x[i] < box.lo.x[i]) dSqr += ((box.lo.x[i]-p.x[i])*(box.lo.x[i]-p.x[i]));
            if (p.x[i] > box.hi.x[i]) dSqr += ((box.hi.x[i]-p.x[i])*(box.hi.x[i]-p.x[i]));
        }
        return dSqr;
    }

    /// Shortest distance between a hyper-box and a point.
    static float_t dist(box_t const& box, point_t const& p) {
        return sqrt(distSqr(box,p) );
    }

    /// Determine whether given 1d-intervals [a0,a1] and [b0,b1] intersect.
    static bool intersects(size_t iDim, double a0, double a1, double b0, double b1) {
        return (std::min(a1, b1) >= std::max(a0, b0) );
    }
};

} /* namespace kdtree */

#endif
