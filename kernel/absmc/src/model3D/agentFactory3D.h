#ifndef ABSMC_AGENT_FACTORY_3D_H
#define ABSMC_AGENT_FACTORY_3D_H

#include <cstddef>
#include <vector>

#include "core/agentFactory.h"
#include "model3D/smc3D.h"
#include "model3D/iel3D.h"
#include "model3D/eel3D.h"
#include "model3D/endo3D.h"
#include "model3D/obstacle3D.h"

namespace absmc {

class AgentFactory3D : public AgentFactory<3> {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<3, double> point_t;

    virtual agent_t* create(AgentTypeId typeId, point_t pos, point_t pos0, double r) const {

    switch (typeId) {
        case tSMC3D:
            return new SMC3D(pos, pos0, r);
        case tIEL3D:
            return new IEL3D(pos, pos0, r);
        case tEEL3D:
            return new EEL3D(pos, pos0, r);
        case tObstacle3D:
            return new Obstacle3D(pos, pos0, r);
        case tEndo3D:
            return new Endo3D(pos, pos0, r);
        default:
            return 0;
    }

    }

    virtual agent_t* create(AgentTypeId typeId, point_t pos, point_t pos0, double r,
                            double drugConc,  double wssOsi, double wssMax, double neoInt, double nitricConc,
                            int age, int clock, int lengthG1,
                            bool baFlag, bool ciFlag, bool selectedForNOFlag) const {

        switch (typeId) {
            case tSMC3D:
                return new SMC3D(pos, pos0, r, drugConc, wssOsi, wssMax, neoInt, nitricConc, age, clock, lengthG1, baFlag, ciFlag, selectedForNOFlag);
            case tIEL3D:
                return new IEL3D(pos, pos0, r);
            case tEEL3D:
                return new EEL3D(pos, pos0, r);
            case tObstacle3D:
                return new Obstacle3D(pos, pos0, r);
            case tEndo3D:
                return new Endo3D(pos, pos0, r, drugConc, wssOsi, wssMax, neoInt, nitricConc, age, clock, lengthG1, baFlag, ciFlag);
            default:
                return 0;
        }

    }

};

} // end namespace absmc

#endif
