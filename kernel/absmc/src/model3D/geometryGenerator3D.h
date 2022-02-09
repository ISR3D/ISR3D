#ifndef ABSMC_GEOMETRY_GENERATOR_3D_H
#define ABSMC_GEOMETRY_GENERATOR_3D_H

#include <vector>

#include "core/geometry.h"

namespace absmc {

class GeometryGenerator3D {
public:
    typedef Point<3, double> point_t;
    typedef Box<3, double> box_t;
    typedef std::vector<size_t> index_vec_t;

    static void generateRegularCuboidDRH(std::vector<point_t> & points, double s, double lX, double lY, double lZ, point_t origin=point_t(0,0,0),
                                         size_t orientation=0, std::vector<index_vec_t> * neighbourhoods=0);

    static void generateSimpleCylinder(std::vector<point_t> & points, double s, double lX, double r1, double r2,
                                       std::vector<index_vec_t> * neighbourhoods=0);
    static void generateSmoothCylinder(std::vector<point_t> & points, double s, double lX, double r1, double r2,
                                       std::vector<index_vec_t> * neighbourhoods=0);

    static void generateSimpleStent(std::vector<point_t> & points, double spacing, double lX, double r, int nStruts, double helixPitch, double strutHeight, double strutWidth);

    static void shufflePositions(std::vector<point_t> & points, double meanShift, double sigmaShift);
    static void shiftPositions(std::vector<point_t> & points, double dX, double dY, double dZ);

    ///s is the spacing between agents; lX, lY and lZ are channel dimensions;
    ///h and w are ridge height and width; start, spacing and count determine the positions of the ridges
    static void generateRidgedChannel(std::vector<point_t> & points, double s, double lX, double lY, double lZ,
                                         double h, double w, double start, double spacing, int count, double wallThickness=0.015,
                                         std::vector<index_vec_t> * neighbourhoods=0);

    ///generates a monolayer with a cons Y coordinate
    static void generateHexagonalMonolayer(std::vector<point_t> & points, double s, double Y,
                                            double lowX, double lowZ, double highX, double highZ,
                                            std::vector<index_vec_t> * neighbourhoods=0);

    static const double EPSILON;

private:

    static void findNbIndices(std::vector<point_t> & points, std::vector<index_vec_t> * neighbourhoods, double range);

    static void completeNbIndicesModuloPeriodicBC(std::vector<point_t> & points, std::vector<index_vec_t> * neighbourhoods, double range,
                                                  size_t iDim, double lo, double hi);
};

} // end namespace absmc

#endif
