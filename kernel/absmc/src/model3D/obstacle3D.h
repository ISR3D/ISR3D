#ifndef ABSMC_OBSTACLE_3D_H
#define ABSMC_OBSTACLE_3D_H

#include <vector>

#include "core/agentBase.h"

namespace absmc {

/// 3d solid obstacle.
class Obstacle3D : public AgentBase<3> {
public:
    typedef AgentBase<3> base_t;
    typedef Point<3, double> point_t;

    Obstacle3D(point_t pos, point_t  pos0, double r)
        : AgentBase<3>(tObstacle3D, pos, pos0, r, point_t(1.0, 1.0, 1.0)), drugConc(0.0)  {
        setAffected(false);
        wssOsi=0.0;
        wssMax=0.0;
    }

    virtual Obstacle3D* clone() const { return new Obstacle3D(*this); }

    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) {
        // CellBase3D::execAgentRule(iAgent, agentsCreated, agentsDeleted);
        // do nothing
    }

    virtual double getScalarByName(std::string const& quantity) const {
        if (quantity=="drugConc") return drugConc;
        else return AgentBase<3>::getScalarByName(quantity);
    }

    double getDrugConc() const { return drugConc; }
    void setDrugConc(double drugConc_) { drugConc = drugConc_; }


    virtual void print(std::ostream & ostr) {
        AgentBase<3>::print(ostr);
        ostr << drugConc << " \t";
        ostr << wssOsi << " \t";
        ostr << wssMax << " \t";
    }

private:
    double drugConc;
};

} // end namespace absmc

#endif
