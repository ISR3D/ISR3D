#ifndef PERCOLATOR_HPP
#define PERCOLATOR_HPP

#include <vector>
#include <stdexcept>
#include <util/typeNames.hpp>

#include "intPoint.hpp"

class Percolator {
private:
    std::vector <int_point> stack;
    std::vector <int_point> seeds;
    int_point domainSize;
    int numSurfaceNodes;

public:
    Percolator() {
    }

    void addSeed(const int_point newSeed) {
        seeds.push_back(newSeed);
    }

    void clearSeeds() {
        seeds.clear();
    }

    std::vector<bool> percolate(const std::vector<int32_t>& nodeTypes, const int_point domainSize_) {
        std::vector<bool> percolated(nodeTypes.size());
        numSurfaceNodes = 0;
        domainSize = domainSize_;

        int x, y, z;

        plantSeeds(percolated, nodeTypes, domainSize);

        while (!stack.empty()) {
            int_point cur = stack.back();
            stack.pop_back();

            x = cur[X];
            y = cur[Y];
            z = cur[Z];

            for (int i = x - 1; i <= x + 1; i++) {
                if (i == -1 || i == domainSize[X]) continue;
                for (int j = y - 1; j <= y + 1; j++) {
                    if (j == -1 || j == domainSize[Y]) continue;
                    for (int k = z - 1; k <= z + 1; k++) {
                        if (k == -1 || k == domainSize[Z]) continue;
                        int_point nb(i,j,k);
                        size_t index = getIndex(nb);
                        // Current node will already be percolated
                        if (!percolated[index]) {
                            percolated[index] = true;
                            if (nodeTypes[index] == FLUID) {
                                stack.push_back(nb);
                            } else {
                                numSurfaceNodes++;
                            }
                        }
                    }
                }
            }
        }

        return percolated;
    }

    int getNumberOfSurfaceNodes() {
        return numSurfaceNodes;
    }

private:
    void plantSeeds(std::vector<bool>& percolated, const std::vector<int32_t>& domain, const int_point domainSize) {
        for (size_t i = 0; i < seeds.size(); ++i) {
            size_t index = getIndex(seeds[i]);

            if (domain[index] != FLUID) {
                throw std::runtime_error("The seed of the percolation filter appears to be solid");
            }
            percolated[index] = true;
            stack.push_back(seeds[i]);
        }
    }

    size_t getIndex(const int_point point) {
        size_t x = point[X];
        size_t y = point[Y];
        size_t z = point[Z];
        size_t index = x*domainSize[1]*domainSize[2] + y*domainSize[2] + z;
        return index;
    }

};

#endif
