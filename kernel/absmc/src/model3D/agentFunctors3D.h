#ifndef ABSMC_AGENT_FUNCTORS_3D_H
#define ABSMC_AGENT_FUNCTORS_3D_H

#include <util/agentTypeId.h>
#include "model3D/smc3D.h"

namespace absmc {

/// Functor setting the agent biologically inactive if it is of type tSMC3D.
class FSetSMC3DInactive  : public std::binary_function<size_t, AgentBase<3>, void> {
public:
    typedef SMC3D agent_t;
    typedef AgentBase<3> base_t;

    void operator()(size_t id, base_t * agent) {
        if (agent->getTypeId() == tSMC3D) {
            SMC3D* smc = dynamic_cast<SMC3D*>(agent);
            smc->setBaFlag(false);
        }
    }
};

/// Functor setting the endothelium to denuded over tSMC3D type cells.
class FDenudeSMC3DEndothelium  : public std::binary_function<size_t, AgentBase<3>, void> {
public:
    typedef SMC3D agent_t;
    typedef AgentBase<3> base_t;

    void operator()(size_t id, base_t * agent) {
        if (agent->getTypeId() == tSMC3D) {
            SMC3D* smc = dynamic_cast<SMC3D*>(agent);
            smc->setSelectedFlag(false);
        }
    }
};

} // end namespace absmc

#endif
