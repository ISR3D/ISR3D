#ifndef ABSMC_INTEGRATOR_H
#define ABSMC_INTEGRATOR_H

#include <vector>

#include "core/integratorController.h"

namespace absmc {

/// Base class for ode integrators.
template<class Agent>
class Integrator {
public:
    typedef Agent agent_t;

    virtual ~Integrator() { }

    /// Integrate the overdamped equation of motion for the system of agents in agentContainer.
    /// The stopping criterion is supplied by the IntegratorController object.
    virtual void integrate(std::vector<agent_t*> & agents, IntegratorController & controller) = 0;
};


} // end namespace absmc

#endif
