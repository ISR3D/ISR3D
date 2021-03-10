#ifndef ABSMC_EEL_RULES_3D_H
#define ABSMC_EEL_RULES_3D_H

#include <vector>
#include <cassert>
#include <iostream>

#include <util/agentTypeId.h>
#include "core/agentRule.h"
#include "model3D/eel3D.h"
#include "core/geometry.h"


namespace absmc {

class EELBreakingRule3D : public AgentRule<EEL3D> {
public:
    typedef EEL3D agent_t;
    typedef AgentBase<3> agent_base_t;
    typedef Point<3, double> point_t;

    EELBreakingRule3D(double maxStrain_, double maxStress_) : maxStrain(maxStrain_), maxStress(maxStress_) { }

    virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

        // if stress is larger than threshold, delete agent
        if (norm<3>(agent.getStress() ) > maxStress) {
            std::cout << "deleting EEL agent: stress=" << norm<3>(agent.getStress()) << " > maxStress" << std::endl;
            agentsDeleted.push_back(&agent);
            return;
        }

        // strain computation ...
        const point_t pos = agent.getPos();
        const point_t pos0 = agent.getPos0();

        // for all neighbours of type EEL3D
        EEL3D::nb_vec_t const& neighbours = agent.getNeighbours();
        const size_t nNb = neighbours.size();
        for (size_t iNb=0; iNb<nNb; iNb++) {

            if (neighbours[iNb]->getTypeId() != tEEL3D) continue;

            // compute initial and current distance and strain
            const point_t nbPos  = neighbours[iNb]->getPos();
            const point_t nbPos0 = neighbours[iNb]->getPos0();

            point_t dist, dist0;
            for (size_t iDim=0; iDim<3; iDim++) {
                dist[iDim]  = pos[iDim]  - nbPos[iDim];
                dist0[iDim] = pos0[iDim] - nbPos0[iDim];
            }
            const double distNorm  = norm<3>(dist);
            const double distNorm0 = norm<3>(dist0);
            const double strain = (distNorm-distNorm0) / distNorm0;
            agent.setStrain(strain); //I have added this line

            // if strain is positive and larger than threshold, delete agent
            if (distNorm0 > 0 && strain > maxStrain) {
                std::cout << "deleting EEL agent: strain=" << strain << " > maxStrain" << std::endl;
                agentsDeleted.push_back(&agent);
                return;
            }
        }
    }

private:
    double maxStrain, maxStress;
};


} // end namespace absmc

#endif
