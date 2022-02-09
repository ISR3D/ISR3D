#ifndef ABSMC_IEL_RULES_3D_H
#define ABSMC_IEL_RULES_3D_H

#include <vector>
#include <cassert>
#include <iostream>

#include <util/agentTypeId.h>
#include "core/agentRule.h"
#include "model3D/iel3D.h"
#include "core/geometry.h"

namespace absmc {

class IELBreakingRule3D : public AgentRule<IEL3D> {
public:
    typedef IEL3D agent_t;
    typedef AgentBase<3> agent_base_t;
    typedef Point<3, double> point_t;

    IELBreakingRule3D(double maxStrain_, double maxStress_) : maxStrain(maxStrain_), maxStress(maxStress_) { }

    virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

        /// if stress is larger than threshold, delete agent
        if (norm<3>(agent.getStress() ) > maxStress) {
            /*not going to kill IEL, but going to transform it into an SMC, so IEL is really a thin layer*/
            point_t const& curPos = agent.getPos();
            point_t const& curPos0 = agent.getPos0();

            agent_base_t* agentB = new SMC3D(curPos, curPos0, agent.getR());
            agentsCreated.push_back(agentB);

            //std::cout << "deleting IEL agent: stress=" << norm<3>(agent.getStress()) << " > maxStress" << std::endl;
            agentsDeleted.push_back(&agent);
            return;
        }

        // strain computation ...
        const point_t pos = agent.getPos();

        // for all bonds of type IEL3D
        IEL3D::bond_vec_t const& bonds = agent.getBonds();
        const size_t nNb = bonds.size();
        size_t num_strains = 0;
        double strain_acc = 0.; //strain accumulator
        for (size_t iNb=0; iNb<nNb; iNb++) {

            if (bonds[iNb].other->getTypeId() != tIEL3D) continue;
            else num_strains += 1;

            // compute initial and current distance and strain
            const point_t nbPos  = bonds[iNb].other->getPos();
            const double distNorm  = norm<3>(pos - nbPos);
            const double strain = (distNorm - bonds[iNb].length0) / bonds[iNb].length0;

            strain_acc += strain;

            // if strain is positive and larger than threshold, break the bond
            if (strain > maxStrain) {
                //std::cout << "breaking IEL bond: pairwise strain=" << strain << " > maxStrain" << std::endl;
                agent.breakBond(bonds[iNb].other);
//                agentsDeleted.push_back(&agent);
//                return;
            }
        }
        if(num_strains > 0)
            agent.setStrain(strain_acc / num_strains);
    }

private:
    double maxStrain, maxStress;
};

} // end namespace absmc

#endif
