#ifndef ABSMC_FORCE_PROTO_H
#define ABSMC_FORCE_PROTO_H

#include "core/geometry.h"
#include "core/euclideanMetrics.h"

namespace absmc {

template<class Agent>
class UnaryForce {
public:
    typedef Point<Agent::nDim, double> point_t;
    typedef EuclideanMetrics<Agent::nDim, double> metrics_t;

    virtual point_t calculateForce(Agent & agent) = 0;

    virtual void accumulate(Agent & agent) {
        agent.addForce(calculateForce(agent));
    }

};

template<class Agent>
class BinaryForce {
public:
    typedef Point<Agent::nDim, double> point_t;

    /// accumulate functions may modify both agents at the same time.
    /// If they do, it's their job to ensure they only do it once each cycle,
    /// since the integrator is likely to pass both (agent1, agent2) and (agent2, agent1) to them.

    /// returns the force vector acting on agent 1
    virtual point_t calculateForce(Agent & agent1, Agent & agent2) = 0;

    //Force > 0 means attraction
    virtual double calculateForce(const double& firstR, const double& secondR, const double& distance) = 0;

    virtual void accumulate(Agent & agent1, Agent & agent2) {
        agent1.addForce(calculateForce(agent1, agent2));
    }

};


} // end namespace absmc

#endif
