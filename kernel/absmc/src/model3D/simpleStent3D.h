#ifndef ABSMC_SIMPLE_STENT_3D_H
#define ABSMC_SIMPLE_STENT_3D_H

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
class SimpleStent3D {
public:
    typedef Point<3, double> point_t;

    /// Create stent from a list of obstacles stored in a text file.
    SimpleStent3D(std::string const& fileName, double cY, double cZ) : cY(cY), cZ(cZ), scalingValue(0.0) {
        // create obstacle agents only for records which match the tObstacle3D typeId
        AgentFileReader<3>::readFile(fileName, AgentFactory3D(), obstacles, tObstacle3D);
    }
    /// Create stent from a vector of obstacles.
    SimpleStent3D(std::vector<Obstacle3D*> newObstacles, double cY, double cZ) : cY(cY), cZ(cZ), scalingValue(0.0) {
        std::copy(newObstacles.begin(), newObstacles.end(), std::back_inserter(obstacles) );
    }

    /// Destructor; does not delete obstacles.
    ~SimpleStent3D () { }

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

    // Places the stent so that it's centered around the central axis f(x)=(x, cY, cZ)
    // and so that its middle point has the x coordinate vesselX/2
    // Afterward, the stent is expanded so that its outer agents are placed at vesselR from the central axis
    void affix(double vesselR, double vesselX, bool scaleR = true) {
        const size_t nObst = obstacles.size();
        point_t center(vesselX * 0.5, cY, cZ);
        point_t min, max;
        findBoundingBox(min, max);

        double maxR = (max[2] - min[2]) / 2; //assuming axisymmetrical stent
        const double ratio = scaleR ? vesselR / maxR - 1.0 : 0.0;
        scalingValue = ratio;

        for (size_t iObst=0; iObst<nObst; iObst++) {
            point_t pos = obstacles[iObst]->getPos();
            //offset so that stent center is in vessel center
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

    double getScalingValue() { return scalingValue; }

private:
    std::vector<AgentBase<3>*> obstacles;
    double cY;
    double cZ;
    double scalingValue;
};


} // end namespace absmc

#endif
