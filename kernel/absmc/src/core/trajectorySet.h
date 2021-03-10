#ifndef ABSMC_TRAJECTORY_SET_H
#define ABSMC_TRAJECTORY_SET_H

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "kdtree/vectorAccessor.h"
#include "kdtree/kdtree.h"
#include "kdtree/kdtree.hh"
#include "core/pointSet.h"
#include "kdtree/euclideanMetrics.h"


namespace absmc {

    using kdtree::Point;
    using kdtree::Box;

    template<size_t nDim>
    class TrajectorySet {
    public:

        typedef AgentBase<nDim> agent_t;
        typedef Point<nDim, double> point_t;
        typedef std::vector<size_t> idx_vec_t;
        typedef EuclideanMetrics<nDim, double> metrics_t;

        TrajectorySet() {};
        ~TrajectorySet() {
            for(auto& step : steps)
                delete step;
        };

        /// Adds a new step. Limitation: steps should be added sequentially
        void addStep(double time, std::vector<point_t>& curCoords) {
            timePoints.push_back(time);
            PointSet<nDim>* temp = new PointSet<nDim>(curCoords);
            steps.push_back(temp);
        }

        void readFromFile(std::string const& fileName) {}

//        getTrajByCurrentPosition


        point_t getTrajByInitPosition(point_t pos0, size_t step, size_t k = 16) {
            //if (k < 2) k = 2;
            idx_vec_t nbs = initCoords.getMClosestIdx(pos0, k);
            std::vector<point_t> nearest (k);
            //put the correct difference in positions into "nearest"
            for (size_t idx = 0; idx < k; idx++) {
                if ((step == 0) || step >= steps.size()){
                    nearest[idx] = point_t(0.,0.,0.);
                    continue;
                }
                nearest[idx] = (*steps[step])[nbs[idx]] - (*steps[step-1])[nbs[idx]];
            }

            // do the interpolation (for k > 1)
            if (k > 1) {
                // calculate distances to the reference points
                std::vector<double> distances (k);
                double distSum = 0.;
                for (size_t idx = 0; idx < k; idx++) {
                    distances[idx] = metrics.dist(initCoords[nbs[idx]], pos0);
                    distSum += distances[idx];
                }

                // calculate weights
                std::vector<double> weights (k);
                double weightsSum = 0.;
                for (size_t idx = 0; idx < k; idx++) {
                    weights[idx] = (1 - distances[idx] / distSum);
                    weightsSum += weights[idx];
                }

                point_t weightedAverage;
                for (size_t idx = 0; idx < k; idx++) {
                    weightedAverage += (weights[idx] / weightsSum) * nearest[idx];
                }
                return weightedAverage;
            }

            return nearest[0];
//            if (step == 0) return (*steps[0])[nbs[0]];
//            if (step >= steps.size()) return point_t(0.,0.,0.);
//            return (*steps[step])[nbs[0]] - (*steps[step-1])[nbs[0]];

        }

        void print() const {
            initCoords.print();

//            for(auto step : steps) {
//                step->print();
//                std::cout << std::endl;
//            }
        }

        size_t stepsCount() const {
            return steps.size();
        }

        void printCount() const {
            std::cout << "Initial coordinates number: " << initCoords.count() << std::endl;
            std::cout << "Steps number: " << steps.size() << std::endl;
//            for(auto step : steps) {
//                std::cout << "Points in step: " << step->count() << std::endl;
//            }
        }

    private:
        std::vector< PointSet<nDim>* > steps;
        std::vector<double> timePoints;
        PointSet<nDim> initCoords;
        metrics_t metrics;
    };


    template<>
    void TrajectorySet<3>::readFromFile (std::string const& fileName) {
        std::ifstream ifstr(fileName.c_str());
        if (!ifstr) {
            std::cout << "Cannot open file " << fileName << " for reading." << std::endl;
            return;
        }
        double curTime = 0.; /// time is not in the current file format, so generate arbitrary sequential timesteps
        std::string line;
        std::vector<double> curValues;
        std::vector<point_t> curPoints;
        std::getline(ifstr, line); /// skip header
        std::getline(ifstr, line); /// read initial positions
        util::tokenizeAndConvert<double>(line, curValues, '\t');
        for(size_t idx = 0; idx <= curValues.size() - 4; idx+=4) {
            curPoints.push_back(point_t(curValues[idx+1], curValues[idx+2], curValues[idx+3]));
        }
        initCoords.initializeFromVector(curPoints);

        std::getline(ifstr, line); /// skip 2nd header
        while (ifstr) {
            curValues.clear();
            curPoints.clear();
            std::getline(ifstr, line); /// read initial positions

            /// an empty line indicates end of file
            if (line.empty() ) break;

            util::tokenizeAndConvert<double>(line, curValues, '\t');
            for(size_t idx = 0; idx <= curValues.size() - 4; idx+=4) {
                curPoints.push_back(point_t(curValues[idx+1], curValues[idx+2], curValues[idx+3]));
            }

            addStep(curTime, curPoints);
            curTime += 1.;
        }
    }

} // end namespace absmc

#endif
