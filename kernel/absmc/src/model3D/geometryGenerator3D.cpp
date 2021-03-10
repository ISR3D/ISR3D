#include <cassert>
#include <iostream>
#include <algorithm>

#include "model3D/geometryGenerator3D.h"
#include <util/random.h>
#include "core/mathUtil.h"
#include <kdtree/vectorAccessor.h>
#include <kdtree/kdtree.h>
#include <kdtree/kdtree.hh>

namespace absmc {

const double GeometryGenerator3D::EPSILON = 1e-4;

//origin is the lowest value for each coordinate
//orientation is the axis normal to the hexagonal planes
void GeometryGenerator3D::generateRegularCuboidDRH(std::vector<point_t> & points, double s, double lX, double lY, double lZ,
                                                    point_t origin, size_t orientation,
                                                    std::vector<index_vec_t> * neighbourhoods)
{
    assert(s>0 && lX>0 && lY>0 && lZ>0);
    if (orientation > 2)
        orientation = 0;
    // spacings
    //or=0
    const double dZ = s;
    const double dY = s * 0.5*sqrt(3.0);
    const double dX = s * sqrt(2.0/3.0);
    //3-2-1

    //or=1
    //const double dZ = s;
    //const double dX = s * 0.5*sqrt(3.0);
    //const double dY = s * sqrt(2.0/3.0);
    //3-1-2

    // or=2
    //const double dX = s;
    //const double dY = s * 0.5*sqrt(3.0);
    //const double dZ = s * sqrt(2.0/3.0);
    //1-2-3

    // number of agents in each direction
    // (neglects agent dimension, first and last agent centre on boundaries)
    // (domain size approximate)
    int nI, nJ, nK; //for orientations 1 and 2 the cuboid is rotated when the points are placed
    //ensuring correct dimensions before rotation
    if (orientation == 0) {
        nI = 1 + (int)floor(lX/dX);
        if (nI % 2 == 0) { nI=nI-1; }   // ensure an odd number of layers
        nJ = 1 + (int)round(lY/dY);
        if (nJ % 2 == 0) { nJ=nJ-1; }   // ensure an odd number of rows
        nK = 1 + (int)round(lZ/dZ);
        if (nK % 2 == 0) { nK=nK-1; }   // ensure an odd number of columns
    }
    if (orientation == 1) {
        nI = 1 + (int)floor(lY/dX);
        if (nI % 2 == 0) { nI=nI-1; }   // ensure an odd number of layers
        nJ = 1 + (int)round(lX/dY);
        if (nJ % 2 == 0) { nJ=nJ-1; }   // ensure an odd number of rows
        nK = 1 + (int)round(lZ/dZ);
        if (nK % 2 == 0) { nK=nK-1; }   // ensure an odd number of columns
    }
    if (orientation == 2) {
        nI = 1 + (int)floor(lZ/dX);
        if (nI % 2 == 0) { nI=nI-1; }   // ensure an odd number of layers
        nJ = 1 + (int)round(lY/dY);
        if (nJ % 2 == 0) { nJ=nJ-1; }   // ensure an odd number of rows
        nK = 1 + (int)round(lX/dZ);
        if (nK % 2 == 0) { nK=nK-1; }   // ensure an odd number of columns
    }

    int n = 0;
    double x, y, z;

    for (int i=0; i<nI; i++) {

        x = i*dX;
        const int layerParity = (i % 2);
        const double y0 = (layerParity == 0) ? 0.0 : 0.5*dZ/sqrt(3.0);

        for (int j=0; j<nJ; j++) {

            y = y0 + j*dY;
            const double z0 = (j % 2 == layerParity) ? 0.0 : 0.5*dZ;

            for (int k=0; k<nK; k++) {
                n++;
                z = z0 + k*dZ;
                if (orientation == 0)
                    points.push_back(point_t(x + origin[0], y + origin[1], z + origin[2]) );
                if (orientation == 1)
                    points.push_back(point_t(y + origin[0], x + origin[1], z + origin[2]) );
                if (orientation == 2)
                    points.push_back(point_t(z + origin[0], y + origin[1], x + origin[2]) );
            }

        }
    }

    // detect all neighbours with equilibrium range + epsilon
    if (neighbourhoods) { findNbIndices(points, neighbourhoods, s+EPSILON); }
}

void GeometryGenerator3D::generateSimpleCylinder(std::vector<point_t> & points, double s, double lX, double r1, double r2,
                                                 std::vector<index_vec_t> * neighbourhoods)
{

    assert(s>0 && lX>0);
    assert(r1>=0 && r2>=0);
    assert(r2 > r1);

    // dimensions of a cuboid containing the cylindrical shell
    const double lY = 2*r2;
    const double lZ = lY;

    // computations as for the drh cuboid
    const double dZ = s;
    const double dY = s * 0.5*sqrt(3.0);

    int nJ = 1 + (int)round(lY/dY);
    if (nJ % 2 == 0) { nJ=nJ-1; }   // ensure an odd number of rows
    int nK = 1 + (int)round(lZ/dZ);
    if (nK % 2 == 0) { nK=nK-1; }   // ensure an odd number of columns

    std::cout << "R1 = " << r1 << std::endl;
    std::cout << "R2 = " << r2 << std::endl;
    std::cout << "lX = " << lX << std::endl;
    std::cout << "lY = " << lY << std::endl;
    std::cout << "lZ = " << lZ << std::endl;

    // coordinates of central axis
    const double cY = 0.5*(nJ-1)*dY;
    const double cZ = 0.5*(nK-1)*dZ;
    std::cout << "cY = " << cY << std::endl;
    std::cout << "cZ = " << cZ << std::endl;

    // generate a regular drh cuboid...
    std::vector<point_t> pointsCuboid;
    std::vector<index_vec_t> neighbourhoodsCuboid;
    if (neighbourhoods) { generateRegularCuboidDRH(pointsCuboid, s, lX, lY, lZ, point_t(0,0,0), 0, &neighbourhoodsCuboid); }
    else { generateRegularCuboidDRH(pointsCuboid, s, lX, lY, lZ, point_t(0,0,0), 0, 0); }
    // ... and shift it such that it is centered around the axis
    shiftPositions(pointsCuboid, 0.0, -cY, -cZ);

    // then copy only points within the cylindrical shell between radii r1 and r2
    const double r1Sqr = r1*r1;
    const double r2Sqr = r2*r2;

    size_t nYP=0, nYN=0, nYZero=0, nZP=0, nZN=0, nZZero=0;

    const size_t nPointsCuboid = pointsCuboid.size();
    for (size_t iPoint=0; iPoint<nPointsCuboid; iPoint++) {

        point_t const& p = pointsCuboid[iPoint];
        const double rSqr = p[1]*p[1]  + p[2]*p[2];

        if (r1Sqr<=rSqr && rSqr<=r2Sqr) {

            // copy point and corresponding neighbourhood if criterion met
            points.push_back(point_t(p) );
            if (neighbourhoods) { neighbourhoods->push_back(neighbourhoodsCuboid[iPoint] );}

            // statistics
            if (p[1] < 0) nYN++;
            if (p[1] > 0) nYP++;
            if (p[1] == 0) nYZero++;
            if (p[2] < 0) nZN++;
            if (p[2] > 0) nZP++;
            if (p[2] == 0) nZZero++;
        }
    }
    std::cout << "n(y<0) = " << nYN << std::endl;
    std::cout << "n(y>0) = " << nYP << std::endl;
    std::cout << "n(y=0) = " << nYZero << std::endl;

    std::cout << "n(z<0) = " << nZN << std::endl;
    std::cout << "n(z>0) = " << nZP << std::endl;
    std::cout << "n(z=0) = " << nZZero << std::endl;
}

void GeometryGenerator3D::generateSmoothCylinder(std::vector<point_t> & points, double s, double lX, double r1, double r2,
                                                 std::vector<index_vec_t> * neighbourhoods)
{
    assert(s>0 && lX>0);
    assert(r1>=0 && r2>=0);
    assert(r2 > r1);

    // approximate dimensions of a cuboid to be transformed into the cylindrical shell;
    // exact dimensions must be computed as (nJ-1)*dY, (nK-1)*dZ.
    const double lY = (r2-r1);
    const double lZ = 2*math::pi * 0.5*(r1+r2);

    // computations as for the drh cuboid

    const double dZ = s;
    const double dX = s * 0.5*sqrt(3.0);
    const double dY = s * sqrt(2.0/3.0);

    int nI = 1 + (int)floor(lX/dX);
    if (nI % 2 == 0) { nI=nI-1; }   // ensure an odd number of layers
    int nJ = 1 + (int)round(lY/dY);
    if (nJ % 2 == 0) { nJ=nJ-1; }   // ensure an odd number of rows
    int nK = 1 + (int)round(lZ/dZ);
    if (nK % 2 == 0) { nK=nK-1; }   // ensure an odd number of columns

    std::cout << "nI=" << nI << endl;
    std::cout << "nJ=" << nJ << endl;
    std::cout << "nK=" << nK << endl;

    std::cout << "lX=" << lX << std::endl;
    std::cout << "(nJ-1)*dY=" << (nJ-1)*dY << std::endl;
    std::cout << "(nK-1)*dZ=" << (nK-1)*dZ << std::endl;

    // generate a regular drh cuboid;
    // code replication from the method above, except that agents in plane z=0 are not generated.
    int n = 0;
    double x, y, z;

    for (int j=0; j<nJ; j++) {

        y = j*dY;
        const int layerParity = (j % 2);
        const double x0 = (layerParity == 0) ? 0.0 : 0.5*dZ/sqrt(3.0);

        for (int i=0; i<nI; i++) {

            x = x0 + i*dX;
            const int maxK = (i % 2 == layerParity) ? nK : nK-1;
            const double z0 = (i % 2 == layerParity) ? 0.0 : 0.5*dZ;

            for (int k=0; k<maxK; k++) {
                n++;
                z = z0 + k*dZ;
                if (z > 0.0) { points.push_back(point_t(x, y, z) ); }
            }

        }
    }

    if (neighbourhoods) {
        // detect all neighbours with equilibrium range + epsilon
        findNbIndices(points, neighbourhoods, s+EPSILON);
        // complete neighbourhoods with neighbours "modulo periodic boundary condition"
        completeNbIndicesModuloPeriodicBC(points, neighbourhoods, s+EPSILON, 2, 0.0, (nK-1)*dZ);
    }

    // transform cuboid into a cylindrical shell by wrapping it around the x-axis.

    // exact middle and inner radius of the wrapping transformation;
    // these may differ from 0.5*(r1+r1) resp. r1 because exact domain dimensions can differ from lY,lZ.
    const double rMean = (nK-1)*dZ / (2.0*math::pi);
    const double rOffs = rMean - 0.5*(nJ-1)*dY;

    std::cout << "R1 = " << r1 << std::endl;
    std::cout << "R2 = " << r2 << std::endl;
    std::cout << "rMean = " << rMean << std::endl;
    std::cout << "rOffs = " << rOffs << std::endl;

    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints; iPoint++) {

        point_t & p = points[iPoint];

        const double phi = p[2] / rMean;
        const double yPrime = (p[1]+rOffs) * cos(phi);
        const double zPrime = (p[1]+rOffs) * sin(phi);

        p[1] = yPrime;
        p[2] = zPrime;
    }

}

void GeometryGenerator3D::generateRidgedChannel (std::vector<point_t> & points, double s, double lX, double lY, double lZ,
                                         double h, double w, double start, double spacing, int count, double wallThickness,
                                         std::vector<index_vec_t> * neighbourhoods)
{
    ///put up the outer walls so that (X,0,0) is the center
    //down
    generateRegularCuboidDRH(points, s, lX, 2 * wallThickness, lZ + 2 * wallThickness,
                            point_t(0, -lY/2 - 2 * wallThickness, -lZ/2 - wallThickness),
                            1, neighbourhoods);
    //up
    generateRegularCuboidDRH(points, s, lX, wallThickness, lZ + 2 * wallThickness,
                            point_t(0, lY/2, - lZ/2 - wallThickness),
                            1, neighbourhoods);
    //right
    generateRegularCuboidDRH(points, s, lX, lY, wallThickness,
                            point_t(0, -lY/2, -lZ/2 - wallThickness),
                            1, neighbourhoods);
    //left
    generateRegularCuboidDRH(points, s, lX, lY, wallThickness,
                            point_t(0, -lY/2, lZ/2),
                            1, neighbourhoods);
    //add ridges
    double ridgeX = start;
    double ridgeThickness = 0.025; //2 * wallThickness;
    for(int i = 0; i < count; i++) {
        //front
        generateRegularCuboidDRH(points, s, ridgeThickness, h - ridgeThickness, lZ,
                            point_t(ridgeX, -lY/2, -lZ/2),
                            0, neighbourhoods);
        //top
        generateRegularCuboidDRH(points, s, w, ridgeThickness, lZ,
                            point_t(ridgeX, -lY/2 + h - ridgeThickness, -lZ/2),
                            1, neighbourhoods);
        //back
        generateRegularCuboidDRH(points, s, ridgeThickness, h - ridgeThickness, lZ,
                            point_t(ridgeX + w - ridgeThickness, -lY/2, -lZ/2),
                            0, neighbourhoods);

        ridgeX = ridgeX + spacing + w;
    }
}

void GeometryGenerator3D::generateHexagonalMonolayer(std::vector<point_t> & points, double s, double Y,
                                            double lowX, double lowZ, double highX, double highZ,
                                            std::vector<index_vec_t> * neighbourhoods)
{
    // spacings
    const double dZ = s;
    const double dX = s * 0.5*sqrt(3.0);

    const double lX = highX - lowX;
    const double lZ = highZ - lowZ;

    // number of agents in each direction
    // (neglects agent dimension, first and last agent centre on boundaries)
    // (domain size approximate)
    int nI = 1 + (int)floor(lX/dX);
    if (nI % 2 == 0) { nI=nI-1; }   // ensure an odd number of layers
    int nK = 1 + (int)round(lZ/dZ);
    if (nK % 2 == 0) { nK=nK-1; }   // ensure an odd number of columns

    double x, y, z;
    y = Y;

    for (int i=0; i<nI; i++) {

        x = lowX + i * dX;
        const double z0 = lowZ + (i % 2 ? 0.0 : 0.5*dZ);

        for (int k=0; k<nK; k++) {
            z = z0 + k * dZ;
            points.push_back(point_t(x, y, z) );
        }

    }

    // detect all neighbours with equilibrium range + epsilon
    if (neighbourhoods) { findNbIndices(points, neighbourhoods, s+EPSILON); }


}


//Generate stent from two sets of spirals, similar to WALLSTENT design
void GeometryGenerator3D::generateSimpleStent(std::vector<point_t> & points, double dr, double lX, double r, int nStruts, double spiralStep)
{
    assert(dr>0 && lX>0);
    assert(r>=0);
    assert(nStruts>0);

    // computations as for the drh cuboid
    const double dX = dr * sqrt(2.0/3.0);
    int nI = 1 + (int)floor(lX/dX);
    if (nI % 2 == 0) { nI=nI-1; }
    const double lXX = nI*dX;

    // number of obstacles in longitudinal (x-) direction
    //const int nII = 1+(int)floor(lXX/s);
    const int nII = 1+2*(int)floor(lXX/dr);

    const double dPhi = 2*math::pi / nStruts;

    // generate 2 * nStruts spiral struts
    double x, y, z;
    for (int iStrut=0; iStrut<nStruts; iStrut++) {
        for (int i=-2; i<3; i++) {
            for (int j=0; j<4; j++) {
                // generate spiral struts at angle phi
                for (int iXX=0; iXX<nII; iXX++) {
                    x = iXX*dr/2.;
                    const double phi1 = (iStrut+0.5) * dPhi + i*dr + (x/spiralStep)* 2*math::pi;
                    const double phi2 = (iStrut+0.5) * dPhi + i*dr - (x/spiralStep)* 2*math::pi;
                    y = (r-j*dr/2) * sin(phi1);
                    z = (r-j*dr/2) * cos(phi1);
                    points.push_back(point_t(x, y, z) );
                    y = (r-j*dr/2) * sin(phi2);
                    z = (r-j*dr/2) * cos(phi2);
                    points.push_back(point_t(x, y, z) );
                }
            }
        }
    }
}

void GeometryGenerator3D::shufflePositions(std::vector<point_t> & points, double meanShift, double sigmaShift)
{
    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints; iPoint++) {

        const double shift = util::Random::getNextGauss(meanShift, sigmaShift);
        const double theta = util::Random::getNext(math::pi);
        const double phi   = util::Random::getNext(2.0*math::pi);

        points[iPoint][0] += shift * sin(theta) * cos(phi);
        points[iPoint][1] += shift * sin(theta) * sin(phi);
        points[iPoint][2] += shift * cos(theta);
    }
}

void GeometryGenerator3D::shiftPositions(std::vector<point_t> & points, double dX, double dY, double dZ)
{
    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints; iPoint++) {

        points[iPoint][0] += dX;
        points[iPoint][1] += dY;
        points[iPoint][2] += dZ;
    }
}

void GeometryGenerator3D::findNbIndices(std::vector<point_t> & points, std::vector<index_vec_t> * neighbourhoods, double range)
{
    assert(neighbourhoods);
    neighbourhoods->clear();

    typedef kdtree::VectorAccessor<3, double, std::vector<point_t> > accessor_t;
    typedef kdtree::KDTree<3, double, accessor_t > tree_t;

    accessor_t acc(points);
    tree_t kdtree(acc);

    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints ; iPoint++) {

        // find indices of all points within a ball of radius range ...
        std::vector<size_t> nbIndices;
        kdtree.findWithinBall(iPoint, range, nbIndices);
        // ... and append the resulting vector of indices to the vector of neighbourhoods
        neighbourhoods->push_back(nbIndices);
    }
}

void GeometryGenerator3D::completeNbIndicesModuloPeriodicBC(std::vector<point_t> & points, std::vector<index_vec_t> * neighbourhoods, double range,
                                                            size_t iDim, double lo, double hi)
{
    assert(neighbourhoods);

    typedef kdtree::VectorAccessor<3, double, std::vector<point_t> > accessor_t;
    typedef kdtree::KDTree<3, double, accessor_t > tree_t;

    // compute rectangular boundary zones of width range at the lower resp. upper boundary
    const double maxVal = std::numeric_limits<double>::max();
    assert(0<=iDim && iDim<=2);

    point_t loBZ1 = point_t(-maxVal, -maxVal, -maxVal);
    loBZ1[iDim] = lo;
    point_t loBZ2 = point_t( maxVal,  maxVal,  maxVal);
    loBZ2[iDim] = lo+range+EPSILON;
    box_t loBZ(loBZ1, loBZ2);

    point_t hiBZ1 = point_t(-maxVal, -maxVal, -maxVal);
    hiBZ1[iDim] = hi-range-EPSILON;
    point_t hiBZ2 = point_t( maxVal,  maxVal,  maxVal);
    hiBZ2[iDim] = hi;
    box_t hiBZ(hiBZ1, hiBZ2);

    // find indices of all agents within boundary zones
    accessor_t acc(points);
    tree_t kdtree(acc);

    std::vector<size_t> idxLoBZ, idxHiBZ;

    kdtree.findWithinBox(loBZ, idxLoBZ);
    kdtree.findWithinBox(hiBZ, idxHiBZ);
    std::cout << idxLoBZ.size() << " points in loBZ\n";
    std::cout << idxHiBZ.size() << " points in hiBZ\n";

    // for all points within the "lower" boundary zone, create a clone "modulo periodic bc",
    // detect neighbours of that clone and add them to the neighbourhood of the original
    const size_t nLoBz = idxLoBZ.size();
    for (size_t iLoBZ=0; iLoBZ<nLoBz; iLoBZ++) {
        // create clone
        const size_t idxPoints = idxLoBZ[iLoBZ];
        point_t const& original = points[idxPoints];
        point_t clone(original);
        clone[iDim] += (hi-lo);
        // detect neighbours of clone
        std::vector<size_t> nbIndicesNew;
        kdtree.findWithinBall(clone, range, nbIndicesNew);
        // add neighbours to neighbourhood of original...
        std::vector<size_t> & nbIndicesAll = (*neighbourhoods)[idxPoints];
        std::copy(nbIndicesNew.begin(), nbIndicesNew.end(), std::back_inserter(nbIndicesAll) );
        // ... and make sure there are no duplicates
        std::sort(nbIndicesAll.begin(), nbIndicesAll.end() );
        nbIndicesAll.erase(std::unique(nbIndicesAll.begin(), nbIndicesAll.end()), nbIndicesAll.end() );
    }

    // same procedure for the "upper" boundary zone
    const size_t nHiBz = idxHiBZ.size();
    for (size_t iHiBZ=0; iHiBZ<nHiBz; iHiBZ++) {
        // create cHine
        const size_t idxPoints = idxHiBZ[iHiBZ];
        point_t const& original = points[idxPoints];
        point_t clone(original);
        clone[iDim] -= (hi-lo);
        // detect neighbours of clone
        std::vector<size_t> nbIndicesNew;
        kdtree.findWithinBall(clone, range, nbIndicesNew);
        // add neighbours to neighbourhood of original...
        std::vector<size_t> & nbIndicesAll = (*neighbourhoods)[idxPoints];
        std::copy(nbIndicesNew.begin(), nbIndicesNew.end(), std::back_inserter(nbIndicesAll) );
        // ... and make sure there are no duplicates
        std::sort(nbIndicesAll.begin(), nbIndicesAll.end() );
        nbIndicesAll.erase(std::unique(nbIndicesAll.begin(), nbIndicesAll.end()), nbIndicesAll.end() );
    }
}

} // end namespace absmc
