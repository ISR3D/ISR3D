#ifndef ABSMC_EEL_3D_H
#define ABSMC_EEL_3D_H

#include <vector>

#include "model3D/cellBase3D.h"
#include "core/agentRule.h"

namespace absmc {

/// 3d EEL agent.
class EEL3D : public CellBase3D {
public:
    typedef AgentBase<3> base_t;
    typedef AgentRule<EEL3D> rule_t;
    typedef Point<3, double> point_t;

    EEL3D(point_t pos, point_t  pos0, double r)
        : CellBase3D(tEEL3D, pos, pos0, r), rule(0) { }

    virtual EEL3D* clone() const { return new EEL3D(*this); }

    void setAgentRule(rule_t *const rule_) {  rule = rule_; }
    rule_t *const getAgentRule() { return rule; }

    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) {
        CellBase3D::execAgentRule(iAgent, agentsCreated, agentsDeleted);
        if (rule) { rule->apply(iAgent, *this, agentsCreated, agentsDeleted); }
    }

private:
    rule_t* rule;
};

} // end namespace absmc

#endif
