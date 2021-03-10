#ifndef ABSMC_SIMPLE_BALLOON_3D_H
#define ABSMC_SIMPLE_BALLOON_3D_H

#include <vector>
#include <algorithm>
#include <cfloat>
#include <cstdlib>

#include "core/agentFileReader.h"
#include "core/agentContainer.h"
#include "core/agentBase.h"
#include "model3D/obstacle3D.h"
#include "model3D/agentFactory3D.h"

namespace absmc {

/// Convenience class to handle all solid obstacles forming a 3d stent.
/// Deployment is modeled by simple radial expansion.
/// Obstacles are not deleted during destruction of a stent object.
/// assuming that the radial direction is in the y-z plane,
/// and the symmetry axis is located at (cY, cZ).
class SimpleBalloon3D {
public:
    typedef Point<3, double> point_t;

    /// Create balloon from a list of obstacles stored in a text file.
    SimpleBalloon3D(std::string const& fileName, double cY, double cZ) : cY(cY), cZ(cZ) {
        // create obstacle agents only for records which match the tObstacle3D typeId
        AgentFileReader<3>::readFile(fileName, AgentFactory3D(), obstacles, tObstacle3D);
    }

    /// Create balloon inside a given arterial geometry
    ///
    SimpleBalloon3D(AgentContainer<3> & agentContainer, double cY, double cZ, double xMin, double xMax, double elementR) :
        cY(cY), cZ(cZ), xMin(xMin), xMax(xMax), elementR(elementR) {

        /// detect smallest free r for each x-value
        const double dx = elementR;
        auto agents = agentContainer.getAgentVector();
        const size_t numX = (xMax - xMin) / dx;
        std::vector<double> freeR;
        freeR.resize(numX+1);

        for(auto const& agent: agents) {
            point_t pos = agent->getPos();
            if (pos[0] < xMin || pos[0] > xMax) continue;
            const size_t xRow = (pos[0] - xMin) / dx;

            // distance from symmetry axis
            const double rSqr = (pos[1]-cY)*(pos[1]-cY) + (pos[2]-cZ)*(pos[2]-cZ);
            const double r = sqrt(rSqr);

            if (r < freeR[xRow] || freeR[xRow] == 0.0) {
                freeR[xRow] = r;
            }
        }
        /// generate balloon agents in a ring shape for each x-value
        /// and add each value to the index
        auto it = std::max_element(std::begin(freeR), std::end(freeR));
        const double angleStep = dx / (*it) ;
        for (size_t xRow = 0; xRow < freeR.size(); xRow++) {
            for (double alpha = 0; alpha < 2 * 3.14159; alpha += angleStep) {
                point_t pos(xMin + xRow * dx, sin(alpha) * (freeR[xRow] - elementR), cos(alpha) * (freeR[xRow] - elementR) );
                Obstacle3D* balloonElem = new Obstacle3D(pos, pos, elementR);
                obstacles.push_back(balloonElem);
            }
        }
    }

    /// Destructor; does not delete obstacles.
    ~SimpleBalloon3D () { }

    /// Return number of obstacles.
    size_t count() const { return obstacles.size(); }

    /// Add all obstacles as agents to the given agent container.
    void exportObstacles(AgentContainer<3> & agentContainer) const {
        agentContainer.add(obstacles);
    }

    /// Increase radial coordinate of all obstacles by dr,

    void expand(double dr) {
        const size_t nObst = obstacles.size();
        for (size_t iObst=0; iObst<nObst; iObst++) {

            point_t pos = obstacles[iObst]->getPos();
            // distance from symmetry axis
            const double rSqr = (pos[1]-cY)*(pos[1]-cY) + (pos[2]-cZ)*(pos[2]-cZ);
            const double r = sqrt(rSqr);

            // obstacles very close to the central axis are possibly not handled correctly;
            // however these should not occur.
            if (r==0) continue;

            // dy/y = dr/r
            pos[1] += (pos[1]-cY) * dr/r;
            // dz/z = dr/r
            pos[2] += (pos[2]-cZ) * dr/r;

            obstacles[iObst]->setPos(pos);
        }
    }

    ///move for dr along the axis X=0, Y=1, or Z=2
    void move(size_t axis, double dr) {
        if (axis > 2) abort();

        const size_t nObst = obstacles.size();
        for (size_t iObst=0; iObst<nObst; iObst++) {
           point_t pos = obstacles[iObst]->getPos();
           pos[axis] += dr;
           obstacles[iObst]->setPos(pos);
        }
    }

    void scale(double ratio, double vesselX) {
        const size_t nObst = obstacles.size();
        point_t min(DBL_MAX,DBL_MAX,DBL_MAX);
        point_t max(-DBL_MAX,-DBL_MAX,-DBL_MAX);
        point_t center(vesselX * 0.5, cY, cZ);

        for (size_t iObst=0; iObst<nObst; iObst++) {
            point_t pos = obstacles[iObst]->getPos();
            for (size_t i = 0; i < 3; ++i) {
                if (pos[i] > max[i]) max[i] = pos[i];
                if (pos[i] < min[i]) min[i] = pos[i];
            }
        }

        for (size_t iObst=0; iObst<nObst; iObst++) {
            point_t pos = obstacles[iObst]->getPos();
            //offset so that balloon center is in vessel center
            for (size_t i = 0; i < 3; ++i) {
                pos[i] += (center[i] - (max[i] + min[i]) * 0.5);
            }
            //r *= vessel R / max R
            //dy/y = dr/r
            pos[1] += (pos[1] - cY) * ratio;
            // dz/z = dr/r
            pos[2] += (pos[2] - cZ) * ratio;

            obstacles[iObst]->setPos(pos);
        }

    }

    void findBoundingBox(point_t& min, point_t& max) {
        min = point_t (DBL_MAX,DBL_MAX,DBL_MAX);
        max = point_t (-DBL_MAX,-DBL_MAX,-DBL_MAX);

        const size_t nObst = obstacles.size();
        for (size_t iObst=0; iObst<nObst; iObst++) {
            point_t pos = obstacles[iObst]->getPos();
            for (size_t i = 0; i < 3; ++i) {
                if (pos[i] > max[i]) max[i] = pos[i];
                if (pos[i] < min[i]) min[i] = pos[i];
            }
        }
    }

    void moveToCenterX (double vesselX) {
        const size_t nObst = obstacles.size();
        point_t center(vesselX * 0.5, cY, cZ);
        point_t min, max;
        findBoundingBox(min, max);

        for (size_t iObst=0; iObst<nObst; iObst++) {
            point_t pos = obstacles[iObst]->getPos();
            //offset so that stent center is in vessel center
            pos[0] += (center[0] - (max[0] + min[0]) * 0.5);
            obstacles[iObst]->setPos(pos);
        }
    }

    void removeBalloon(AgentContainer<3> & agentContainer) {
        std::set< AgentBase<3>* > setToRemove(obstacles.begin(), obstacles.end());
        agentContainer.remove(setToRemove);
    }

private:
    std::vector<AgentBase<3>*> obstacles;
    double cY;
    double cZ;
    double xMin, xMax, elementR;
};


} // end namespace absmc

#endif
