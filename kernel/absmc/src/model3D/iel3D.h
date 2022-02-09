#ifndef ABSMC_IEL_3D_H
#define ABSMC_IEL_3D_H

#include <vector>

#include "model3D/cellBase3D.h"
#include "core/agentRule.h"

namespace absmc {

/// 3d IEL agent.
class IEL3D : public CellBase3D {
public:
    typedef AgentBase<3> base_t;
    typedef AgentRule<IEL3D> rule_t;
    typedef Point<3, double> point_t;

    IEL3D(point_t pos, point_t  pos0, double r)
        : CellBase3D(tIEL3D, pos, pos0, r), rule(0) { }

    virtual IEL3D* clone() const { return new IEL3D(*this); }

    void setAgentRule(std::shared_ptr<rule_t> const rule_) {  rule = rule_; }
    std::shared_ptr<rule_t> const getAgentRule() { return rule; }

    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) {
        CellBase3D::execAgentRule(iAgent, agentsCreated, agentsDeleted);
        if (rule) { rule->apply(iAgent, *this, agentsCreated, agentsDeleted); }
    }

private:
    std::shared_ptr<rule_t> rule;
};

} // end namespace absmc

#endif
