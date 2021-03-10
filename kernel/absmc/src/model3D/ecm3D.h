#ifndef ABSMC_ECM_3D_H
#define ABSMC_ECM_3D_H

#include <vector>

#include "model3D/cellBase3D.h"
#include "core/agentRule.h"

namespace absmc {

/// 3d IEL agent.
class ECM3D : public CellBase3D {
public:
    typedef AgentBase<3> base_t;
    typedef AgentRule<ECM3D> rule_t;
    typedef Point<3, double> point_t;

    //enum CycleState { G0, G1, SG2M };

    ECM3D(point_t pos, point_t  pos0, double r)
        : CellBase3D(tECM3D, pos, pos0, r),
          rule(0)  { }

    virtual ECM3D* clone() const { return new ECM3D(*this); }

    void setAgentRule(rule_t *const rule_) {  rule = rule_; }
    rule_t *const getAgentRule() { return rule; }

    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) {
        CellBase3D::execAgentRule(iAgent, agentsCreated, agentsDeleted);
        if (rule) { rule->apply(iAgent, *this, agentsCreated, agentsDeleted); }
    }
    virtual double getScalarByName(std::string const& quantity) const {
        //if (quantity=="state") return static_cast<int>(state);
        //else if (quantity=="age") return age;
        return CellBase3D::getScalarByName(quantity);
    }

    /// Increase (or decrease) the cell's notional radius by the given growth factor.
    void grow(double growthFactor) { setR( getR()*growthFactor ); }

    virtual void print(std::ostream & ostr) {
        CellBase3D::print(ostr);
        //ostr << ciFlag << " \t";
    }

private:
    rule_t* rule;

};

} // end namespace absmc

#endif
