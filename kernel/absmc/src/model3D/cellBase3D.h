#ifndef ABSMC_CELL_BASE_3D_H
#define ABSMC_CELL_BASE_3D_H

#include "core/agentBase.h"
#include "core/agentRule.h"

namespace absmc {

/// Base class for cells of the 3d model. Defines variables for quantities computed by other models of the COAST ISR CxA.
/// Mobility is set to (1.0, 1.0, 1.0) for all 3d cells.
class CellBase3D : public AgentBase<3> {
public:
    typedef AgentBase<3> base_t;
    typedef AgentRule<CellBase3D> rule_t;
    typedef Point<3, double> point_t;

    /// Constructor.
    CellBase3D(AgentTypeId typeId, point_t pos, point_t  pos0, double r)
        : AgentBase<3>(typeId, pos, pos0, r, point_t(1.0, 1.0, 1.0)), drugConc(0.0), neoDist(0.0), nitricConc(0.0), rule(0) {

        wssOsi=0.0;
        wssMax=0.0;
    }

    CellBase3D(AgentTypeId typeId, point_t pos, point_t  pos0, double r, double drugConc,  double wssOsi_, double wssMax_, double neoDist, double nitricConc)
    : AgentBase<3>(typeId, pos, pos0, r, point_t(1.0, 1.0, 1.0)), drugConc(drugConc), neoDist(neoDist), nitricConc(nitricConc), rule(0) {

        wssOsi = wssOsi_;
        wssMax = wssMax_;
    }

    virtual double getScalarByName(std::string const& quantity) const {
        if (quantity=="drugConc") return drugConc;
        else if (quantity=="neoDist") return neoDist;
        else return AgentBase<3>::getScalarByName(quantity);
    }

    void setNeoDist(double neoDist_) { neoDist = neoDist_; }
    double getNeoDist() const { return neoDist; }

    double getDrugConc() const { return drugConc; }
    void setDrugConc(double drugConc_) { drugConc = drugConc_; }

    double getNitricConc() const { return nitricConc; }
    void setNitricConc(double nitricConc_) { nitricConc = nitricConc_; }

    virtual void print(std::ostream & ostr) {
        AgentBase<3>::print(ostr);
        ostr << drugConc << " \t";
        ostr << wssOsi << " \t";
        ostr << wssMax << " \t";
        ostr << neoDist << " \t";
    }

    void setAgentRule(rule_t *const rule_) {  rule = rule_; }
    rule_t *const getAgentRule() { return rule; }

    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) {
        if (rule) { rule->apply(iAgent, *this, agentsCreated, agentsDeleted); }
    }

private:
    double drugConc;
    double neoDist;
    double nitricConc;   /// Nitric Oxide Concentration in Pa
    rule_t* rule;
};

} // end namespace absmc

#endif
