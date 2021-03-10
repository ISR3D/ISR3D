#ifndef ABSMC_SMC_3D_H
#define ABSMC_SMC_3D_H

#include <vector>

#include "model3D/cellBase3D.h"
#include "core/agentRule.h"
#include <util/random.h>

namespace absmc {

/// 3d smooth muscle cell.
class SMC3D : public CellBase3D {
public:
    typedef AgentBase<3> base_t;
    typedef AgentRule<SMC3D> rule_t;
    typedef Point<3, double> point_t;

    enum CycleState { G0, G1, SG2M };
    SMC3D(point_t pos, point_t  pos0, double r)
    : CellBase3D(tSMC3D, pos, pos0, r), rule(0), state(G0), age(0), clock(0), baFlag(true), ciFlag(false), selectedForNOFlag(true),
          lengthG1(static_cast<int>(util::Random::getNextGauss(SMC3D::lengthG1Mean, SMC3D::lengthG1StdDev)) ) { }
    // end of debug

    SMC3D(point_t pos, point_t  pos0, double r, double drugConc,  double wssOsi, double wssMax, double neoDist, double nitricConc,
            int age, int clock, int lengthG1,
            bool baFlag, bool ciFlag, bool selectedForNOFlag)
    : CellBase3D(tSMC3D, pos, pos0, r, drugConc, wssOsi, wssMax, neoDist, nitricConc), rule(0), state(G0),
    age(age), clock(clock), baFlag(baFlag), ciFlag(ciFlag), selectedForNOFlag(selectedForNOFlag),
    lengthG1(lengthG1) { }

    virtual SMC3D* clone() const { return new SMC3D(*this); }

    void setAgentRule(rule_t *const rule_) {  rule = rule_; }
    rule_t *const getAgentRule() { return rule; }

    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) {
        CellBase3D::execAgentRule(iAgent, agentsCreated, agentsDeleted);
        if (rule) { rule->apply(iAgent, *this, agentsCreated, agentsDeleted); }
    }

    virtual double getScalarByName(std::string const& quantity) const {
        if (quantity=="state") return static_cast<int>(state);
        else if (quantity=="age") return age;
        else if (quantity=="clock") return clock;
        else if (quantity=="ciFlag") return ciFlag;
        else if (quantity=="lengthG1") return lengthG1;
        else if (quantity=="baFlag") return baFlag;
        else if (quantity=="noFlag") return noFlag;
        else if (quantity=="selectedForNOFlag") return selectedForNOFlag;
        else return CellBase3D::getScalarByName(quantity);
    }

    CycleState getState() const { return state; }
    void setState(CycleState state_) { state = state_; }

    int getAge() const { return age; }
    void setAge(int age_) { age = age_; }

    int getClock() const { return clock; }
    void setClock(int clock_) { clock = clock_; }

    bool getBaFlag() const { return baFlag;}
    void setBaFlag(bool baFlag_) { baFlag = baFlag_; }

    bool getCiFlag() const { return ciFlag; }
    void setCiFlag(bool ciFlag_) { ciFlag = ciFlag_; }

    bool getNOFlag() const { return noFlag; }
    void setNOFlag(bool noFlag_) { noFlag = noFlag_; }

    bool getSelectedFlag() const { return selectedForNOFlag; }
    void setSelectedFlag(bool selectedForNOFlag_) { selectedForNOFlag = selectedForNOFlag_; }

    /// Increase the cell's notional radius by the given growth factor.
    void grow(double growthFactor) { setR( getR()*growthFactor ); }

    virtual void print(std::ostream & ostr) {
        CellBase3D::print(ostr);
        ostr << age << " \t";
        ostr << clock << " \t";
        ostr << lengthG1 << " \t";
        ostr << baFlag << " \t";
        ostr << ciFlag << " \t";
        ostr << selectedForNOFlag << " \t";

    }

private:
    rule_t* rule;
    CycleState  state;          // cell cycle state
    int         age;            // number of mitotic cycles completed
    int         clock;          // cell cycle clock
    bool        baFlag;           // biological activity flag
    bool        ciFlag;         // contact inhibition flag
    bool        noFlag;         // Nitric Oxide flag
    bool        selectedForNOFlag; //check if an SMC is elegible for NO flag or not

public:
    static const int lengthG1Mean   = 16;
    static const int lengthG1StdDev = 2 ;  // mean and width of distribution of G1 phase length (in hours)
    static const int lengthSG2M     = 16;  // fixed length of SG2M phase (in hours)
    const int lengthG1;                    // length of G1 phase for this cell

};

} // end namespace absmc

#endif
