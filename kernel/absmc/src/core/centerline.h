#ifndef ABSMC_CENTERLINE_H
#define ABSMC_CENTERLINE_H

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include <kdtree/vectorAccessor.h>
#include <kdtree/kdtree.h>
#include <kdtree/kdtree.hh>
#include "core/pointSet.h"

namespace absmc {

    using kdtree::Point;
    using kdtree::Box;

    template<size_t nDim>
    class Centerline : public PointSet<nDim> {
    public:

        typedef PointSet<nDim> base_;
        typedef Point<nDim, double> point_t;

        Centerline() {};
        ~Centerline() {};

        /// Fill centerline with zeros (Y, Z) for a range of X
        void resetToZero (double minX, double maxX, size_t num) {
            base_::points.resize(num, point_t());
            double curX = minX;
            /// put points at both minX and maxX
            double step = (maxX - minX) / (num - 1);
            for(auto point : base_::points) {
                point[0] = curX;
                curX += step;
            }
            base_::initializeFromVector(base_::points);
        }

        /// Fill centerline with (cY, cZ) for a range of X
        void resetToAxis (double minX, double maxX, double cY, double cZ, size_t num) {
            base_::points.resize(num, point_t());
            double curX = minX;
            /// put points at both minX and maxX
            double step = (maxX - minX) / (num - 1);
            for(auto point : base_::points) {
                point[0] = curX;
                point[1] = cY;
                point[2] = cZ;
                curX += step;
            }
            base_::initializeFromVector(base_::points);
        }

        /// True if a "direction" vector originating from point "coords" goes in the same direction as the vector from "coord" to centerline (< pi/2 angle difference)
        bool isFacingCenterline(point_t coords, point_t direction) const {
            point_t vector_to_centerline = base_::getClosest(coords) - coords;
            double dotproduct = vector_to_centerline * direction;
            return (dotproduct > 0);
        }

    };

} // end namespace absmc

#endif
