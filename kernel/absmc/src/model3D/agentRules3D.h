#ifndef ABSMC_AGENT_RULES_3D_H
#define ABSMC_AGENT_RULES_3D_H

#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <map>

#include "core/agentRule.h"
#include "model3D/cellBase3D.h"
#include "core/geometry.h"
#include "core/euclideanMetrics.h"
#include <util/agentTypeId.h>

namespace absmc {

    /// Agent radius change rule. Used for bending the vessel during deployment.
    class AgentCurveVesselRule3D : public AgentRule< CellBase3D > {
    public:
        typedef CellBase3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;


        AgentCurveVesselRule3D(double curvatureRadius_, const AgentContainer<3> &agents, size_t nSteps_)
        : curvatureRadius(curvatureRadius_), nSteps(nSteps_) {
            nAgents = agents.count();
            targetRadius.resize(nAgents);
            originalRadius.resize(nAgents);
            for (size_t iAgent=0; iAgent < nAgents; iAgent++) {
                double yCoord = agents.getAgentVector()[iAgent]->getScalarByName("posY");
                double agentR = agents.getAgentVector()[iAgent]->getR();
                this->targetRadius[iAgent] = agentR * (curvatureRadius + yCoord) / curvatureRadius;
                this->originalRadius[iAgent] = agentR;
            }
        }

        virtual void apply(size_t iAgent, agent_t & agent,
                            std::vector<agent_base_t*> & agentsCreated,
                            std::vector<agent_base_t*> & agentsDeleted) const
        {
            // don't curve the stent
            if (agent.getTypeId() != tObstacle3D) {
                agent.setR(agent.getR() + (this->targetRadius[iAgent] - this->originalRadius[iAgent]) / nSteps);
            }
        }

    private:
        const double curvatureRadius;
        size_t nAgents;
        const size_t nSteps;
        std::vector<double> targetRadius;
        std::vector<double> originalRadius;
    };

} // end namespace absmc

#endif
