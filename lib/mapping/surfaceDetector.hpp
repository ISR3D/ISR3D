#ifndef SURFACE_DETECTOR_HPP
#define SURFACE_DETECTOR_HPP

#include <vector>
#include "intPoint.hpp"
#include "percolator.hpp"

class SurfaceDetector {

    /**
     * Detect surface nodes
     * @param nodeTypes 3d mesh of node types
     * @return An array of surface nodes [ [i0,j0,k0],[i1,j1,k1],... ]
     */
public:
    virtual std::vector<int_point> detectSurface (const std::vector<int32_t>& nodeTypes, int_point domainSize) = 0;

};

class NoSurface : public SurfaceDetector {
public:
    std::vector<int_point> detectSurface (const std::vector<int32_t>& nodeTypes, int_point domainSize) {
        std::vector<int_point> res(0);
        return res;
    }
};


class PercolationDetector : public SurfaceDetector {
private:
    const int_point seed;
    Percolator percolator;

    /**
     * The seed is a point where the percolation starts.
     * @param seed
     */
public:
    PercolationDetector(int_point seed_) :
    seed(seed_){
        percolator.addSeed(seed);
    }

    std::vector<int_point> detectSurface (const std::vector<int32_t>& nodeTypes, int_point domainSize) {
//        System.out.println("");
//        System.out.println("[SURFACE]: Seed=" + seed[X] + ", " + seed[Y] + ", " + seed[Z]);
//        System.out.println("[SURFACE]: SeedType=" + domain[seed[X]][seed[Y]][seed[Z]]);
//        System.out.println("");

        //Get surface
        std::vector<bool> percolated = percolator.percolate(nodeTypes, domainSize);

        // Form a vector of surface nodes
        std::vector<int_point> surfaceNodes;
        for (int i = 0; i < domainSize[X]; i++) {
            for (int j = 0; j < domainSize[Y]; j++) {
                for (int k = 0; k < domainSize[Z]; k++) {
                    size_t index = i*domainSize[1]*domainSize[2] + j*domainSize[2] + k;
                    if (nodeTypes[index] == SOLID && percolated[index]) {
                        surfaceNodes.push_back(int_point (i,j,k));
                    }
                }
            }
        }

//        System.out.println("[SURFACE]: Surface nodes = " + surfaceNodes.size());
        return surfaceNodes;
    }

};

#endif
