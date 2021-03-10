#ifndef SPHERE_VOXELIZER_HPP
#define SPHERE_VOXELIZER_HPP


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <list>
#include <stdexcept>
#include <util/logger.h>
#include <util/typeNames.hpp>

using util::logger;

#include "morphologyTools3D.hpp"
#include "sphere.hpp"
#include "surfaceDetector.hpp"
#include "intPoint.hpp"
#include "percolator.hpp"

class SphereVoxelizer {
private:
    typedef Point<3, double> point_t;
    typedef EuclideanMetrics<3, double> metrics_t;

    const std::string outPath;
    std::vector<Sphere> shapes;
    // std::vector<std::vector<std::vector<unsigned char> > > nodeTypes;
    // std::vector<std::vector<std::vector<int> > > shapeReference;
    std::vector<int32_t> nodeTypes;
    std::vector< std::vector<size_t>* > shapeReference;
    std::vector<int_point> surfaceLocation;
    const double dx;
    const std::vector<double> domain;
    const double smoothRadius;
    SurfaceDetector& surfaceDetector;
    bool outerSideIsPercolated;

    void discretize();
    void detectSurface();
    void mapToShapes(double value, int i, int j, int k);
    std::vector<double> getAndResetShapeMapping();
    std::vector<double> getAndResetShapeMappingMin();
    std::vector<double> getAndResetShapeMappingMax();
    void clearShapeReference();

    void checkVoxelization() {
        if (nodeTypes.empty()) {
            throw std::runtime_error ("Please, call the voxelize method before.");
        }
    }
    void smooth() {
        int smoothNumber = static_cast<int>(ceil(smoothRadius / dx)); //TODO: check if it is ceil or floor.
        nodeTypes = mathematicalMorphologyTools3D::open(smoothNumber, nodeTypes, sizeX, sizeY, sizeZ, SOLID);
    }

public:
    const size_t sizeX;
    const size_t sizeY;
    const size_t sizeZ;

    static NoSurface noSurface;
    SphereVoxelizer(double dx, std::vector<double> domain, double smoothRadius = 0.0, SurfaceDetector& surfaceDetector = noSurface);
    ~SphereVoxelizer () {
        dispose();
    }
    std::vector<Sphere> coordinates2Shapes(std::vector<double> coords);
    //void addToShapes(std::vector<double> coords);
    bool nodeExists(int i, int j, int k) const;
    void addShapeReference(int i, int j, int k, size_t idx);
    void percolateOutside();
    void setMaskedDomain(std::vector<int32_t> domain);
    void writeSurfaceToFile(std::string filename);
    std::vector<double> getLumenArea() const;
    std::vector<double> mapSurfaceAverage(std::vector<double> data, std::vector<int32_t> mask);
    std::vector<double> mapDomain(std::vector<double> values, double threshold);
    std::vector<double> mapDomainMax(std::vector<double> values, double threshold);

    static int coord2Idx(double coord, double origin, double dx) {
        return static_cast<int>(ceil((coord - origin) / dx));
    }
    static std::vector<int> coord2Idx(std::vector<double> coords, std::vector<double> domain, double dx) {
        std::vector<int> toReturn;

        toReturn.push_back(static_cast<int>(round((coords[0] - domain[0])/dx)));
        toReturn.push_back(static_cast<int>(round((coords[1] - domain[1])/dx)));
        toReturn.push_back(static_cast<int>(round((coords[2] - domain[2])/dx)));
        toReturn.push_back(static_cast<int>(round((coords[3] - domain[0])/dx)));
        toReturn.push_back(static_cast<int>(round((coords[4] - domain[1])/dx)));
        toReturn.push_back(static_cast<int>(round((coords[5] - domain[2])/dx)));

        return toReturn;
    }

    static point_t idx2Coord(int_point idx, std::vector<double> domain, double dx) {
        point_t coords;
        for (int i = 0; i < 3; i++) {
            coords[i] = domain[i] + idx[i] * dx + 0.5 * dx;
        }
        return coords;
    }

    static point_t idx2Coord(int i, int j, int k, std::vector<double> domain, double dx) {
        int_point idx (i, j, k);
        return idx2Coord(idx, domain, dx);
    }

    /*
     * Return the size of the nodeTypes vector.
     */
    int_point getSize() {
        return int_point(sizeX, sizeY, sizeZ);
    }
    /**
     * Voxelize a previously read list of cells.
     * The <tt>dispose()</tt> method is automatically called.
     * @param coords An 1d array of cells parameters: [x0,y0,z0,r0,x1,y1,z1,r1, ...}
     */
    void voxelize() {
        discretize();
        if (smoothRadius > 0.0) {
            smooth();
        }
        detectSurface();
    }
    void voxelizeNoSurface() {
        discretize();
        if (smoothRadius > 0.0) {
            smooth();
        }
    }
    std::vector< std::vector<size_t>* > getShapeReference() {
        return this->shapeReference;
    }
    std::vector<int_point> getSurface() {
        return this->surfaceLocation;
    }
    void setSurface(std::vector<int_point> surface) {
        this->surfaceLocation = surface;
    }
    std::vector<Sphere> getShapes() {
        return this->shapes;
    }
    void setShapes(std::vector<double> coords) {
        this->shapes = coordinates2Shapes(coords);
    }

    void setShapes(std::vector<Sphere>& shapes) {
        this->shapes = shapes;
    }

    void setShapes(std::list<Sphere>& shapes) {
        this->shapes = { std::begin(shapes), std::end(shapes) };
    }

    void setShapeReference(std::vector< std::vector<size_t>* > shapeRef) {
        this->shapeReference = shapeRef;
    }

    /// Clears the vectors used by the voxelizer
    void dispose() {
       shapes.clear();
       nodeTypes.clear();
       surfaceLocation.clear();
       clearShapeReference();
    }
    /**
     * Returns the voxelized, smoothed domain by reference.
     * @return A 3d array of integers indicating the voxel type. //TODO: refer to proper class
     */
    std::vector<int32_t> getSmoothDomain() {
        checkVoxelization();
        return nodeTypes;
    }
};

/**
* Create a new sphere voxelizer instance. All the parameters are in physical units.
* @param dx The mesh delta x
* @param domain The two opposite corner of the domain bounding box: [x0, y0, z0, x1, y1, z1]
* @param smoothRadius The smoothness radius. Set it to zero to avoid smoothing,
* @param surfaceDetector The class used to detect surface
*/
SphereVoxelizer::SphereVoxelizer(double dx, std::vector<double> domain, double smoothRadius, SurfaceDetector& surfaceDetector) :
    dx(dx), domain(domain), smoothRadius(smoothRadius), surfaceDetector(surfaceDetector),
    // abs in case of negative values
    sizeX (abs(coord2Idx(domain[3], domain[0], dx))),
    sizeY (abs(coord2Idx(domain[4], domain[1], dx))),
    sizeZ (abs(coord2Idx(domain[5], domain[2], dx))) {

    outerSideIsPercolated = false;
    shapeReference.resize(sizeX * sizeY * sizeZ);
    // TODO: path has to be passed specifically
    //outPath = instance.get_setting_as<std::string> ("muscle_data_out_dir");
    //logger().setFile(outPath + "SphereVoxelizer.log");
}

/**
*  Converts a list of coordinates into <B>Sphere</B> shapes.
* @param coords in the form of (x,y,z,r).
* @return
*/
std::vector<Sphere> SphereVoxelizer::coordinates2Shapes(std::vector<double> coords) {
   // final int numOfShapes = coords.length / 4;
   const int numOfShapes = coords.size() / 4;

   std::vector<Sphere> s;
   for (int i = 0; i < numOfShapes; i++) {
       double x = coords[i * 4];
       double y = coords[i * 4 + 1];
       double z = coords[i * 4 + 2];
       double r = coords[i * 4 + 3];
        // 42 is not a real type; I think we don't use this function anyway, but let's check
       s.push_back(Sphere(x, y, z, r, 42));
   }
   return s;
}

/**
* Adds a list of coordinates to shapes list.
* Copies the shape array, so use sparingly.
* @param coords in the form of (x,y,z,r).
* @return
*/
//void SphereVoxelizer::addToShapes(std::vector<double> coords) {
//   //const int numOfNewShapes = coords.size() / 4;
//   //int numOfOldShapes = shapes.size();
//   std::vector<Sphere> newShapes = coordinates2Shapes(coords);
//   std::vector<Sphere> s;

//   for (int i : shapes) {
//       s.push_back(i);
//   }
//   for (int j : newShapes) {
//       s.push_back(j);
//   }

//   this->shapes = s;
//
//   ///won't work unless the 2nd vector is shapes
//   shapes.insert(shapes.end(), coords.begin(), coords.end());
//}

bool SphereVoxelizer::nodeExists(int i, int j, int k) const {
    return (i >= 0) && (j >= 0) && (k >= 0) &&
           ((size_t)i < sizeX) && ((size_t)j < sizeY) && ((size_t)k < sizeZ);
}

/**
* Discretize the shapes.
*/
void SphereVoxelizer::discretize() {
    nodeTypes.clear();
    nodeTypes.resize(sizeX*sizeY*sizeZ);
    clearShapeReference();
    shapeReference.resize(sizeX*sizeY*sizeZ);

    for (size_t idx = 0; idx < shapes.size(); ++idx) {
        Sphere& s = shapes[idx];
        std::vector<double> boundBox = s.getBoundingBox();
        std::vector<int> voxelShapeBoundBox = coord2Idx(boundBox, domain, dx);
        for (int i = voxelShapeBoundBox[0]; i <= voxelShapeBoundBox[3]; i++) {
            for (int j = voxelShapeBoundBox[1]; j <= voxelShapeBoundBox[4]; j++) {
                for (int k = voxelShapeBoundBox[2]; k <= voxelShapeBoundBox[5]; k++) {
                    if (nodeExists(i, j, k) && (s.getSize() > 0)) {
                        if (s.isInsideBox(idx2Coord(i, j, k, domain, dx), dx)) {
                            nodeTypes[i*sizeY*sizeZ + j*sizeZ + k] = SOLID;
                        }
                        addShapeReference(i, j, k, idx);
                    }
                }
            }
        }
    }
}

void SphereVoxelizer::clearShapeReference() {
    for (auto& shapePointer : shapeReference) {
        delete shapePointer;
    }
    shapeReference.clear();
}

void SphereVoxelizer::addShapeReference(int i, int j, int k, size_t idx) {
    if (shapeReference[i*sizeY*sizeZ + j*sizeZ + k] == nullptr) {
        shapeReference[i*sizeY*sizeZ + j*sizeZ + k] = new std::vector<size_t>;
    }
    shapeReference[i*sizeY*sizeZ + j*sizeZ + k]->push_back(idx);
}

void SphereVoxelizer::detectSurface() {
   logger() << "SURFACE DETECTION" << std::endl;
   surfaceLocation = surfaceDetector.detectSurface(nodeTypes, getSize());
   logger() << "Found: " << surfaceLocation.size() << " elements" <<std::endl;
   // Add solidMask on the outer side of the vessel
}

void SphereVoxelizer::percolateOutside() {
    logger() << "Percolating the outside of the vessel" << std::endl;
    /// Add solidMask on the outer side of the vessel
    if (!outerSideIsPercolated) {
        Percolator percolator;
        // Add seeds in all 'horizontal' corner rows
        for (size_t x = 0; x < sizeX; x++) {
            percolator.addSeed(int_point(x, 0, 0));
            percolator.addSeed(int_point(x, 0, sizeZ - 1));
            percolator.addSeed(int_point(x, sizeY - 1, 0));
            percolator.addSeed(int_point(x, sizeY - 1, sizeZ - 1));
        }
        std::vector<bool> solidMask = percolator.percolate(nodeTypes, getSize());
        for (size_t i = 0; i < sizeX; i++) {
            for (size_t j = 0; j < sizeY; j++) {
                for (size_t k = 0; k < sizeZ; k++) {
                    size_t index = i*sizeY*sizeZ + j*sizeZ + k;
                    if (solidMask[index] && nodeTypes[index] == FLUID) { //aka percolated from the outside & not solid => outside the vessel
                        nodeTypes[index] = STATIC;
                    }
                }
            }
        }
        outerSideIsPercolated = true;
    }
}

void SphereVoxelizer::setMaskedDomain(std::vector<int32_t> domain) {
    this->nodeTypes = domain;
    for (size_t i = 0; i < sizeX; i++) {
        for (size_t j = 0; j < sizeY; j++) {
            for (size_t k = 0; k < sizeZ; k++) {
                if (domain[i*sizeY*sizeZ + j*sizeZ + k] == STATIC) {
                    outerSideIsPercolated = true;
                    return;
                }
            }
        }
    }
}

std::vector<double> SphereVoxelizer::getLumenArea() const {
    std::vector<double> lumenArea (sizeX);
    for (size_t i = 0; i < sizeX; i++) {
        for (size_t j = 0; j < sizeY; j++) {
            for (size_t k = 0; k < sizeZ; k++) {
                if (this->nodeTypes[i*sizeY*sizeZ + j*sizeZ + k] == FLUID) {
                     lumenArea[i] += 1;
                }
            }
        }
    }
    for(auto& area : lumenArea) {
        area *= dx*dx;
    }
    return lumenArea;
}

void SphereVoxelizer::writeSurfaceToFile(std::string filename) {
    try{
       std::ofstream file;
       file.open(filename);
       for (size_t ijk = 0; ijk < surfaceLocation.size(); ijk++) {
           int i = surfaceLocation[ijk][0];
           int j = surfaceLocation[ijk][1];
           int k = surfaceLocation[ijk][2];

           file << (domain[0] + dx * i) << " " << (domain[1] + dx * j) << " " << (domain[2] + dx * k) << std::endl;
       }
       file.close();
    } catch (std::runtime_error& e) {
        logger() << "RUNTIME ERROR: " << e.what() << std::endl;
        return;
    }
}


   /// various mappers

   /**
    * Maps surface values to the surface cells. Value from the closest cell is mapped.
    * @param values The values to map to the cell [v0,v1...]
    * @param coordinates The corresponding mesh coordinates [ [x0,y0,z0], [x1,y1,z1], ...]
    * @return A value for each cell.
    */
//    public double[] mapSurface(double[] data, char[] mask) {
//        for (int ijk = 0; ijk < surfaceLocation.length; ijk++) {
//            int i = surfaceLocation[ijk][0];
//            int j = surfaceLocation[ijk][1];
//            int k = surfaceLocation[ijk][2];

//            // int index = ii + jj * sizeX + kk * sizeX * sizeY;
//            int indexpalab = sizeZ * (i * sizeY + j) + k;
//            //outbb.write(shearpalabosData[indexpalab]+" "+boolmask[indexpalab]+" ");

//            double closestValue = 0.0;
//            int minDist = Integer.MAX_VALUE;
//            for (int x = -1; x <= 1; x++) {
//                int xAbs = x == -1 ? 1 : x;
//                for (int y = -1; y <=  1; y++) {
//                    int yAbs = y == -1 ? 1 : y;
//                    for (int z = -1; z <= 1; z++) {
//                        if (x == 0 && y == 0 && z == 0) continue;
//                        int zAbs = z == -1 ? 1 : z;
//                         int newindex = indexpalab + sizeZ*(x*sizeY + y) + z;
//                        if (xAbs+yAbs+zAbs < minDist && newindex>=0 && newindex < mask.length && (mask[newindex] == FLUID)) {
//                            closestValue = data[newindex];
//                            minDist = xAbs + yAbs + zAbs;
//                        }
//                    }
//                }
//            }
//            this->mapToShapes(closestValue, i, j, k);
//        }
//        return getAndResetShapeMapping();
//    }

       /**
    * Maps surface values to the surface cells. Maximal value among neighbouring liquid cells is mapped.
    * @param values The values to map to the cell [v0,v1...]
    * @param coordinates The corresponding mesh coordinates [ [x0,y0,z0], [x1,y1,z1], ...]
    * @return A value for each cell.
    */
//    public double[] mapSurfaceMax(double[] data, char[] mask) {
//        for (int ijk = 0; ijk < surfaceLocation.length; ijk++) {
//            int i = surfaceLocation[ijk][0];
//            int j = surfaceLocation[ijk][1];
//            int k = surfaceLocation[ijk][2];

//            int indexpalab = sizeZ * (i * sizeY + j) + k;
//            double maxValue = - Double.MAX_VALUE;
//            for (int x = -1; x <= 1; x++) {
//                for (int y = -1; y <=  1; y++) {
//                    for (int z = -1; z <= 1; z++) {
//                        if (x == 0 && y == 0 && z == 0) continue;
//                         int newindex = indexpalab + sizeZ*(x*sizeY + y) + z;
//                        if (newindex>=0 && newindex < mask.length && (mask[newindex] == FLUID) && data[newindex] > maxValue) {
//                            maxValue = data[newindex];
//                        }
//                    }
//                }
//            }
//            this->mapToShapes(maxValue, i, j, k);
//        }
//        return getAndResetShapeMapping();
//    }

/**
* Maps surface values to the surface cells. Average liquid-cell value in neighbourhood is mapped.
* @param values The values to map to the cell [v0,v1...]
* @param coordinates The corresponding mesh coordinates [ [x0,y0,z0], [x1,y1,z1], ...]
* @return A value for each cell.
*/
std::vector<double> SphereVoxelizer::mapSurfaceAverage(std::vector<double> data, std::vector<int32_t> mask) {
    for (size_t ijk = 0; ijk < surfaceLocation.size(); ijk++) {
        int i = surfaceLocation[ijk][0];
        int j = surfaceLocation[ijk][1];
        int k = surfaceLocation[ijk][2];

        // int index = ii + jj * sizeX + kk * sizeX * sizeY;
        int indexpalab = sizeZ * (i * sizeY + j) + k;

        for (int x = -1; x <= 1; x++) {
            for (int y = -1; y <=  1; y++) {
                for (int z = -1; z <= 1; z++) {
                    if (x == 0 && y == 0 && z == 0) continue;
                    int newindex = indexpalab + sizeZ*(x*sizeY + y) + z;
                    if (newindex>=0 && (size_t)newindex < mask.size() && (mask[newindex] == FLUID)) {
                        this->mapToShapes(data[newindex], i, j, k);
                    }
                }
            }
        }
    }
    return getAndResetShapeMapping();
}


/**
* Map domain values to the cells by average.
* @param values The domain values arranged in a 3d mesh.
* @return A value for each cell.
*/
std::vector<double> SphereVoxelizer::mapDomain(std::vector<double> values, double threshold) {
   //For every node in the domain
   for (size_t i = 0; i < sizeX; i++) {
       for (size_t j = 0; j < sizeY; j++) {
           for (size_t k = 0; k < sizeZ; k++) {
               if (values[i*sizeY*sizeZ + j*sizeZ + k]>threshold) {
                   mapToShapes(values[i*sizeY*sizeZ + j*sizeZ + k], i, j, k);
               }
           }
       }
   }
   return getAndResetShapeMapping();
}

/**
* Map domain values to the cells by maximum.
* @param values The domain values arranged in a 3d mesh.
* @return A value for each cell.
*/
std::vector<double> SphereVoxelizer::mapDomainMax(std::vector<double> values, double threshold) {
   //For every node in the domain
   for (size_t i = 0; i < sizeX; i++) {
       for (size_t j = 0; j < sizeY; j++) {
           for (size_t k = 0; k < sizeZ; k++) {
               if (values[i*sizeY*sizeZ + j*sizeZ + k]>threshold) {
                   mapToShapes(values[i*sizeY*sizeZ + j*sizeZ + k], i, j, k);
               }
           }
       }
   }
   return getAndResetShapeMappingMax();
}

void SphereVoxelizer::mapToShapes(double value, int i, int j, int k) {
   std::vector<size_t>* nodeShapes = shapeReference[i*sizeY*sizeZ + j*sizeZ + k];
   if (nodeShapes != nullptr && !nodeShapes->empty()) {
       for (size_t p = 0; p < nodeShapes->size(); p++) {
           //Add the corresponding Value
           shapes[nodeShapes->at(p)].addMappingValue(value);
       }
   }
}

std::vector<double> SphereVoxelizer::getAndResetShapeMapping() {
   std::vector<double> result(shapes.size());
   //For every cell
   for (size_t i = 0; i < shapes.size(); i++) {
       //Put the average value in the result vector
       result[i] = shapes[i].getMappingAverage();
       //Reset the cell
       shapes[i].resetMapping();
   }
   return result;
}

std::vector<double> SphereVoxelizer::getAndResetShapeMappingMin() {
    std::vector<double> result(shapes.size(), 0.);
    //For every cell
    for (size_t i = 0; i < shapes.size(); i++) {
        //Put the average value in the result vector
        result[i] = shapes[i].getMappingMin();
        //Reset the cell
        shapes[i].resetMapping();
    }
   return result;
}

std::vector<double> SphereVoxelizer::getAndResetShapeMappingMax() {
    std::vector<double> result(shapes.size(), 0.);
    //For every cell
    for (size_t i = 0; i < shapes.size(); i++) {
        //Put the max value in the result vector
        result[i] = shapes[i].getMappingMax();
        //Reset the cell
        shapes[i].resetMapping();
    }
   return result;
}

#endif
