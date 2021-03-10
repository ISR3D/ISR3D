#ifndef ABSMC_POINT_SET_H
#define ABSMC_POINT_SET_H

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "kdtree/vectorAccessor.h"
#include "kdtree/kdtree.h"
#include "kdtree/kdtree.hh"

namespace absmc {

    using kdtree::Point;
    using kdtree::Box;

    template<size_t nDim>
    class PointSet  {
    public:

        typedef AgentBase<nDim> agent_t;
        typedef Point<nDim, double> point_t;
        typedef kdtree::VectorAccessor<nDim, double, std::vector<point_t> > accessor_t;
        typedef kdtree::KDTree<nDim, double, accessor_t > tree_t;
        typedef std::vector<size_t> idx_vec_t;

        PointSet() {}
        ~PointSet() {}

        PointSet(std::vector<point_t> & pointVector) {
            initializeFromVector(pointVector);
        }

        void readFromFile(std::string const& fileName) {
            std::ifstream ifstr(fileName.c_str());
            if (!ifstr) {
                std::cout << "Cannot open file " << fileName << " for reading." << std::endl;
                return;
            }
            point_t pos;
            while (ifstr) {
                for (size_t iDim=0; iDim<nDim; iDim++) ifstr >> pos[iDim];
                points.push_back(pos);
            }
            initializeFromVector(points);
        }

        void writeToFile(std::string const& fileName) {
            std::ofstream ofstr(fileName.c_str());
            if (!ofstr) {
                std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
                return;
            }
            point_t pos;
            for (auto& pos : points) {
                for (size_t iDim=0; iDim<nDim; iDim++) ofstr << pos[iDim] << '\t';
                ofstr << std::endl;
            }
            ofstr.close();
        }

        void initializeFromVector (std::vector<point_t> const& pointVector) {
            points = pointVector;
            accessor_t acc(points);
            delete point_tree;
            point_tree = new tree_t(acc);
        }

        point_t const& getClosest(point_t target) const {
            idx_vec_t temp;
            point_tree->findMNearest(target, 1, temp);
            return points[temp[0]];
        }

        idx_vec_t getMClosestIdx(point_t target, size_t m = 1) const {
            idx_vec_t temp;
            point_tree->findMNearest(target, m, temp);
            return temp;
        }

        point_t const& operator[](size_t i) const {
            return points[i];
        }

        void print() const {
            size_t n = 0;
            for(auto point : points) {
                std::cout << n << ": " << point << std::endl;
            }
        }

        /// Count all points
        size_t count() const {
            return points.size();
        }

    protected:

        /// keep points in a k-d tree for ease of geometric stuff
        std::vector<point_t> points;
        tree_t* point_tree=nullptr;
    };

} // end namespace absmc

#endif
