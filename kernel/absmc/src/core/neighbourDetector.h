#ifndef ABSMC_NEIGHBOUR_DETECTOR_H
#define ABSMC_NEIGHBOUR_DETECTOR_H

#include <cstddef>
#include <vector>
#include <algorithm>

namespace absmc {

/// Base class for neighbour detection modules. Subclasses implement different neighbourhood detection policies.
template<class Agent>
class NeighbourDetector {
public:
    typedef Agent agent_t;

    virtual ~NeighbourDetector() { }

    /// Initialise neighbours for all agents contained in the agent vector.
    virtual void init(std::vector<agent_t*> agents) = 0;

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a slow time scale. Called after deployment step and after application of biological ruleset.
    virtual void updateSlowTimeScale(std::vector<agent_t*> agents) = 0;

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a fast time scale. Called at every iteration of the integrator.
    virtual void updateFastTimeScale(std::vector<agent_t*> agents, double displMax) = 0;

};


} // end namespace absmc

#endif
