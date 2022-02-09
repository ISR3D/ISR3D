#ifndef ABSMC_SIMPLE_BALLOON_3D_H
#define ABSMC_SIMPLE_BALLOON_3D_H

#include <vector>
#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <cmath>

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
    SimpleBalloon3D(double elementR, double length, double balloonR) :
        cY(0.0), cZ(0.0), elementR(elementR) {

        /// detect smallest free r for each x-value
        const double spacing = elementR;
        const double dx = spacing * 0.5*sqrt(3.0);
        const size_t numX = (length) / dx + 1;

        /// generate balloon agents in a ring shape for each x-value
        /// and add each value to the index
        const double angleStep = spacing / balloonR;
        for (size_t xRow = 0; xRow < numX; xRow++) {
            for (double alpha = (xRow % 2 == 0) ? 0. : angleStep / 2; alpha < 2. * math::pi; alpha += angleStep) {
                point_t pos(xRow * dx, sin(alpha) * (balloonR), cos(alpha) * (balloonR) );
                Obstacle3D* balloonElem = new Obstacle3D(pos, pos, elementR);
                obstacles.push_back(balloonElem);
            }
        }

        /// add spherical sections at the ends of the stent
        const double numEdge = (balloonR / dx);
        for (size_t xRow = 0; xRow < numEdge; xRow++) {
            const double edgeR = sqrt((balloonR * balloonR) - (xRow * dx) * (xRow * dx));
            const double angleStepEdge = spacing / edgeR;
            for (double alpha = (xRow % 2 == 0) ? 0. : angleStepEdge / 2; alpha < 2. * math::pi; alpha += angleStepEdge) {
                point_t pos1(0. - xRow * dx, sin(alpha) * (edgeR), cos(alpha) * (edgeR) );
                point_t pos2(length + xRow * dx, sin(alpha) * (edgeR), cos(alpha) * (edgeR) );
                Obstacle3D* balloonElem1 = new Obstacle3D(pos1, pos1, elementR);
                Obstacle3D* balloonElem2 = new Obstacle3D(pos2, pos2, elementR);
                obstacles.push_back(balloonElem1);
                obstacles.push_back(balloonElem2);
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

    void moveToCenterX (double centerX) {
        const size_t nObst = obstacles.size();
        point_t center(centerX, cY, cZ);
        point_t min, max;
        findBoundingBox(min, max);

        for (size_t iObst=0; iObst<nObst; iObst++) {
            point_t pos = obstacles[iObst]->getPos();
            //offset so that stent center is in vessel center
            pos[0] += (center[0] - (max[0] + min[0]) * 0.5);
            obstacles[iObst]->setPos(pos);
        }
    }

    void remove(AgentContainer<3> & agentContainer) {
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
