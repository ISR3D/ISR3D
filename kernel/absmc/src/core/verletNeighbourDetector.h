#ifndef ABSMC_VERLET_NEIGHBOUR_DETECTOR_H
#define ABSMC_VERLET_NEIGHBOUR_DETECTOR_H

#include <cstddef>
#include <cassert>
#include <vector>

#include "core/neighbourDetector.h"
#include "core/neighbourCache.h"
#include "core/neighbourCacheWrapper.h"

namespace absmc {

/// Neighbourhood detection policy effectively equivalent to the "current" policy:
/// compute neighbourhood records from scratch at every call to update() both on fast and slow time scales.
/// However this implementation is more efficient, since a Verlet technique is used to reduce the number of exhaustive searches,
/// i.e. the detection range is larger than the specified cut-off range. Use of this policy makes sense only
/// in conjunction with a force class that cuts off interaction at the specified cut-off range.
template<class Agent>
class VerletNeighbourDetector : public NeighbourDetector<Agent> {
public:
    typedef Agent agent_t;

    VerletNeighbourDetector(double rCutOff_)
        : rCutOff(rCutOff_), rVerlet(2*rCutOff) { }

    VerletNeighbourDetector(double rCutOff_, double rVerlet_)
        : rCutOff(rCutOff_), rVerlet(rVerlet_) { }

    /// Initialise neighbours for all agents contained in the agent vector.
    virtual void init(std::vector<agent_t*> agents) {
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, rVerlet, true);
    }

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a slow time scale. Called after deployment step and after application of biological ruleset.
    virtual void updateSlowTimeScale(std::vector<agent_t*> agents) {
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, rVerlet, true);
    }

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a fast time scale. Called at every iteration of the integrator.
    virtual void updateFastTimeScale(std::vector<agent_t*> agents, double displMax) {
        // Verlet trick: perform exhaustive search only if the maximal displacement of any agent
        // since the last exhaustive update is larger than the thickness of the Verlet skin.
        assert(displMax >= 0.0);
        displSum += displMax;

        if (2*displSum >= rVerlet-rCutOff) {
            NeighbourCache<agent_t> nbCache(agents);
            NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, rVerlet, true);
            displSum = 0.0;
        }
    }

private:
    const double rCutOff;   //  interaction cut-off radius
    const double rVerlet;   //  Verlet skin radius
    double displSum;        // maximal displacement of any agent since the last exhaustive update
};

} // end namespace absmc

#endif
