#ifndef ABSMC_SMC_RULES_3D_H
#define ABSMC_SMC_RULES_3D_H

#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <map>
#include <limits>

#include <util/agentTypeId.h>
#include "core/agentRule.h"
#include "model3D/smc3D.h"
#include "model3D/endo3D.h"
#include "model3D/ecm3D.h"
#include "model3D/endoRules3D.h"
#include "core/geometry.h"
#include "core/euclideanMetrics.h"


namespace absmc {

    /// The contact inhibition rule only sets SMC3D::ciFlag, which is then read by the cell cycle rule.
    /// Before applying the rule to any cell, make sure the neighbour cache pointer is valid.
    class SMCContactInhibitionRule3D : public AgentRule<SMC3D> {
    public:
        typedef SMC3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;

        // the last argument is indeed of type pointer to a pointer... see context and you will understand
        SMCContactInhibitionRule3D(double wSMC_, double wIEL_, double wEEL_, double wObstacle_, double rangeFactor, double contactThreshold_,
                                   NeighbourCache<agent_base_t> ** nbCachePtr_)
        : wSMC(wSMC_), wIEL(wIEL_), wEEL(wEEL_), wObstacle(wObstacle_), rangeFactor(rangeFactor), contactThreshold(contactThreshold_), nbCachePtr(nbCachePtr_) { }

        virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

            NeighbourCache<agent_base_t>* nbCache = *nbCachePtr;
            assert(nbCache);
            std::vector<agent_base_t*> const& agents = nbCache->getAgentVector(); //was this previous non-pointer thing consuming a lot?
            std::vector<size_t> nbIndices;
            double R1 = agent.getR();
            //Arbitrary cutoff; assuming there are no cells bigger than 1.7 R1 nearby;
            //true if the dispersion of initial radii is relatively small
            nbCache->getNbIndices(iAgent, 2.7 * rangeFactor * R1, nbIndices);

            // count contacts of different types
            size_t numSMC = 0, numIEL = 0, numEEL = 0, numObstacle = 0, numECM = 0;

            const size_t nNb = nbIndices.size();
            for (size_t iNb=0; iNb<nNb; iNb++) {
                agent_base_t const* nb = agents[nbIndices[iNb] ];
                //always count the EEL if it's in the neighbourhood
                //this is not really a physical interaction;
                //more like a boundary condition to stop the growth on the outer surface
                if (nb->getTypeId() == tEEL3D) numEEL++;
                //For other types of agents, check if they're close enough to interact
                if (metrics.distSqr(agent.getPos(),nb->getPos()) <= ((R1 + nb->getR())*rangeFactor)*
                                                                    ((R1 + nb->getR())*rangeFactor)) {
                    if (nb->getTypeId() == tSMC3D) numSMC++;
                    if (nb->getTypeId() == tIEL3D) numIEL++;
                    if (nb->getTypeId() == tObstacle3D) numObstacle++;
                    if (nb->getTypeId() == tECM3D) numECM++;
                }
            }

            // actual contact inhibition criterion
            double contactWeight = wSMC*numSMC + wIEL*numIEL + wEEL*numEEL + wObstacle*numObstacle + wSMC*numECM; //TODO: separate weight for ECM
            if (contactWeight > contactThreshold) {
                agent.setCiFlag(true);
            //    return;
            }
            //agent.setCiFlag(false);
        }

    private:
        const double wSMC, wIEL, wEEL, wObstacle, rangeFactor, contactThreshold;
        NeighbourCache<agent_base_t> ** nbCachePtr;
        EuclideanMetrics<3, double> metrics;
    };



    /// For endothelium-covered cells, calculate the NO concentration and set the flag according to WSS
    class SMCNitricOxideRule3D : public AgentRule<SMC3D> {
    public:
        typedef SMC3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;

        SMCNitricOxideRule3D() { }

        virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated,
                            std::vector<agent_base_t*> & agentsDeleted) const {

            const double wss = agent.getWssMax()*10.;// value in dynes/cm2 from value in Pa
            double nitricValue = 0.;
            //check if endothelium grows over this cell
            if (agent.getSelectedFlag()) {
                if (wss >= 0.0 && wss <= 10.5) {
                    nitricValue = 0.000366*wss;    // since the reference fig 2b, Acta physiologica paper, 2009, 197, 99–106
                                                //starts from wss=10.5, so we extracted a equation of a linear line between 0
                                                // and 10.5 to capture the effects of low wss.
                                                //nitricValue in nmol/mm3 aka millimolar
                } else if (wss > 10.5) //agent.wss_max > 1.05
                    nitricValue = (0.00052*(wss*wss))-(0.011*wss)+0.062;    // with reference to fig 2b, above A. Phys. publication
                                                                            //>=0.0038 fro WSS > 10.5 dynes/cm2
            }
            agent.setNitricConc(nitricValue);

            //if (nitricValue >= 0.0004)  //threshold value, equivalent to WSS_Max in agent > 0.109 Pa
                                        //or wss = 1.09 dynes/cm^2
                                        //Hannan's thesis lists threshold value as 0.001 nmol/mm3 = 1000 nM
                                        //the reference paper (Coneski 2012) lists 400 nM as cell cycle arrest value though
                                        //(and 1000 nM as "full cycle arrest", so that must be the upper bound)

            if (nitricValue >= 0.001) //threshold value, equivalent to WSS_Max in agent > 0.273 Pa
                                        //or wss = 2.73 dynes/cm^2
                                        //This goes with Hannan's version
                                        //and is closer to the rough experimental value of 0.5 Pa for arrest
                agent.setNOFlag(true);
        }

    };




    /// HMOX1 rule, based on work by Federico Vozzi et al., in InSilc project
    // works in SMCs directly exposed to the flow (no EC coverage)
    // HMOX expression is linear between 0 and 1 dynes/cm2, as well as between 1 and 20 dynes/cm2
    // Expression at 20 dynes/cm2 corresponds to 100% HMOX activity
    // A percent of SMCs equal to the value of HMOX activity is inhibited (e.g. at 20% of activity 20% of exposed SMCs stop proliferating, and at 100% proliferation stops completely)
    // This rule reuses NOFlag and should be included AFTER the NO rule. If this rule is here to stay, it should get its own flag.
    class HMOX1Rule3D : public AgentRule<SMC3D> {
    public:
        typedef SMC3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;
        const double conc0dyn = 414.7;
        const double conc1dyn = 12383.;
        const double conc20dyn = 35597.;

        //probability range from 0 to 1, per hour
        HMOX1Rule3D(double dt_)
        : dt(dt_)
        { }

        virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated,
                            std::vector<agent_base_t*> & agentsDeleted) const {

            //check if endothelium grows over this cell
            if (!agent.getSelectedFlag()) {
                double wss = agent.getWssMax()*10; // value in dynes/cm2 from value in Pa
                if (wss == 0.) return;

                double HMOX1expression = 0.; //linear between 0 and 1 dynes/cm2, as well as between 1 and 20 dynes/cm2
                if (wss > 0. && wss < 1.)
                    HMOX1expression = conc0dyn + (wss / 1.) * (conc1dyn - conc0dyn);
                if (wss >= 1. && wss < 20.)
                    HMOX1expression = conc1dyn + ((wss - 1.) / (20. - 1.)) * (conc20dyn - conc1dyn);
                if (wss >= 20.)
                    HMOX1expression = conc20dyn;
                double HMOX1activity = HMOX1expression / conc20dyn; //from 0 to 1, linearly proportionate to expression, max at 20dyn

                if (util::Random::getNext(1.0) <= HMOX1activity)
                    agent.setNOFlag(true); // should be a separate HMOX1 flag
            }
        }

    private:
        double dt;
    };




    /// SMC cell cycle rule. The contact inhibition rule must be applied before to set contact inhibition flags.
    class SMCCellCycleRule3D : public AgentRule<SMC3D> {
    public:
        typedef SMC3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;

        SMCCellCycleRule3D(double drugConcThreshold_, double dt_)
        : drugConcThreshold(drugConcThreshold_), dt(dt_)  {
        }


        virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

            const int tCheckPoint = 0;  // G0 checkpoint (at beginning of cell cycle)
            const int tMitosis = agent.lengthG1 + SMC3D::lengthSG2M; //hours
            const int t = agent.getClock(); //hours
            const double growthFactor = pow(2.0, 1.0/(3.0*agent.lengthG1/dt));
            const int intdt = (int)(dt + 0.5);

            // at G0 check point or in G1
            if (t < agent.lengthG1) {
                // check for: biological activity, drug concentration, wss threshold, contact inhibition.
                // if any criterion for prohibition is met, G0.
                if (!agent.getBaFlag() || agent.getDrugConc() >= drugConcThreshold || agent.getCiFlag() || agent.getNOFlag())
                {
                    //return to G0 and stop proliferation
                    if (t > 0) {
                        agent.setState(SMC3D::G0);
                        agent.setClock(tCheckPoint);
                        //agent.setBaFlag(false);
                    }
                }
                // otherwise enter G1 and increase clock
                else {
                    agent.setState(SMC3D::G1);
                    agent.setClock(t+intdt);
                    agent.grow(growthFactor);
                }
            }
            // S/G2/M
            else {
                if (t < tMitosis) { // cell in S/G2/M
                    agent.setState(SMC3D::SG2M);
                    agent.setClock(t+intdt);
                }
                else { // cell at M: perform mitosis

                    // random displacement from position of parent cell
                    const double newR = agent.getR()/ pow(2.0,1.0/3.0);
                    const double shift = 0.4*newR;
                    const double theta = util::Random::getNext(math::pi);
                    const double phi   = util::Random::getNext(2.0*math::pi);

                    const double dX = shift * sin(theta) * cos(phi);
                    const double dY = shift * sin(theta) * sin(phi);
                    const double dZ = shift * cos(theta);

                    point_t const& curPos = agent.getPos();
                    const point_t posA(curPos[0]+dX, curPos[1]+dY, curPos[2]+dZ);
                    const point_t posB(curPos[0]-dX, curPos[1]-dY, curPos[2]-dZ);

                    agent.setPos(posA);
                    agent.setR(newR);
                    agent.setState(SMC3D::G0);
                    agent.setClock(0);
                    agent.setAge(agent.getAge()+1);
                    agent.clearNeighbours();

                    SMC3D* agentB = new SMC3D(posB, posB, newR);
                    agentB->setAgentRule(agent.getAgentRule() );
                    agentB->setSelectedFlag(agent.getSelectedFlag());

                    agentB->setStrain(agent.getStrain() );

                    agent.clearBonds();
                    agent.addBond(agentB);
                    agentB->copyNeighbours(agent);

                    /// add bonds to the surrounding cells, up to 12 bonds total
                    agent.addNeighboursToBonds(true, 12, newR * 2.5);
                    agentB->addNeighboursToBonds(true, 12, newR * 2.5);

                    agentsCreated.push_back(agentB);

                }
            }
        }

    private:
        const double drugConcThreshold;
        const double dt;
    };




    class SMCNecrosisRule3D : public AgentRule<SMC3D> {
    public:
        typedef SMC3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;

        SMCNecrosisRule3D(double maxStress_) : maxStress(maxStress_) { }

        virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

            // if stress is larger than threshold, delete agent
            if (norm<3>(agent.getStress() ) > maxStress) {
                // std::cout << "SMC necrosis rule, deleting agent: stress=" << norm<3>(agent.getStress()) << " > maxStress" << std::endl;
                agentsDeleted.push_back(&agent);
                return;
            }
        }

    private:
        double maxStress;
    };




    class SMCBondBreakingRule3D : public AgentRule<SMC3D> {
    public:
        typedef SMC3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;

        SMCBondBreakingRule3D(double maxStrain_) : maxStrain(maxStrain_) { }

        virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

            // strain computation
            const point_t pos = agent.getPos();

            // for all bonds
            auto& bonds = agent.getBonds();
            for (auto& bond : bonds) {

                // compute equilibrium and current distance and strain
                const point_t nbPos  = bond.other->getPos();
                const double distNorm  = norm<3>(pos - nbPos);
                const double equilibriumDist = bond.other->getR() + agent.getR();
                const double strain = (distNorm-equilibriumDist) / equilibriumDist;

                // if strain is larger than the threshold, break the bond
                if (strain > maxStrain) {
                    //std::cout << "breaking SMC bond: pairwise strain=" << strain << " > maxStrain" << std::endl;
                    agent.breakBond(bond.other);
                }
            }

        }

    private:
        double maxStrain;
    };




    class ECMProductionRule3D : public AgentRule<SMC3D> {
    public:
        typedef SMC3D agent_t;
        typedef AgentBase<3> agent_base_t;
        typedef Point<3, double> point_t;

        ECMProductionRule3D(const double maxStrain_, const double ecmProductionProbability_,  Centerline<3> & centerline_) :
            maxStrain(maxStrain_), ecmProductionProbability(ecmProductionProbability_), centerline(&centerline_) { }

        virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

            if(agent.getStrain() == 0) {
                // strain computation ...
                const point_t pos = agent.getPos();

                // for all bonds of type SMC3D
                SMC3D::bond_vec_t const& bonds = agent.getBonds();
                const size_t nNb = bonds.size();
                size_t num_strains = 0;
                double strain_acc = 0.; //strain accumulator
                for (size_t iNb=0; iNb<nNb; iNb++) {

                    // if (bonds[iNb].other->getTypeId() != tSMC3D) continue;
                    if (bonds[iNb].length0 < 0.001) continue;

                    // compute initial and current distance and strain
                    const point_t nbPos  = bonds[iNb].other->getPos();
                    const double distNorm  = norm<3>(pos - nbPos);
                    const double strain = (distNorm - bonds[iNb].length0) / bonds[iNb].length0;
                    if (strain > 0.) {
                        strain_acc += strain;
                        ++num_strains;
                    }
                }

                double strain = num_strains > 0 ? strain_acc / (double)num_strains : 0.;

                if (std::isinf(strain))
                    strain = 0.0;

                agent.setStrain(strain);
            }

            //produce ECM
            if(agent.getStrain() > maxStrain) {
                //if an inhibition criterion is met, return
                if (!agent.getBaFlag() || agent.getCiFlag() || agent.getNOFlag()) {
                    return;
                }

                //ensure the motility of synthetic cells; clean up the neighbours outside the interaction raduis
                agent.clearNeighbours();

                if (util::Random::getNext(1.0) <= ecmProductionProbability) {
                    //spawn ECM;
                    const double newR = 0.015;
                    const double shift = 0.8 * newR;
                    const point_t pos = agent.getPos();

                    //ECM placement away from X axis (away from the centerline)
                    const point_t center  = centerline->getClosest(pos);
                    const point_t yzProjection = (pos - center);
                    const double yzProjectionNorm  = norm<3>(yzProjection);

                    // place the cell in a semi-spere centered around the line
                    // here, we convert the line to polar coordinates, and add a random component (-pi/2, +pi/2)

                    const double xVec = yzProjection[0] / yzProjectionNorm;
                    const double yVec = yzProjection[1] / yzProjectionNorm;
                    const double zVec = yzProjection[2] / yzProjectionNorm;

                    // "zVec / norm", but the norm = 1
                    const double theta = acos(zVec) - math::pi/2. + util::Random::getNext(math::pi);
                    const double phi   = atan2(yVec, xVec) - math::pi/2. + util::Random::getNext(math::pi);

                    const double dXC = shift * sin(theta) * cos(phi);
                    const double dYC = shift * sin(theta) * sin(phi);
                    const double dZC = shift * cos(theta);

                    const point_t posC(pos[0]+dXC, pos[1]+dYC, pos[2]+dZC);

                    ECM3D* agentC = new ECM3D(posC, posC, newR);

                    agentC->copyNeighbours(agent);
                    /// add bonds to the surrounding cells, up to 12 bonds total
                    agentC->addNeighboursToBonds(true, 12, newR * 2.5);

                    //agentC->setAgentRule(ecmRule);
                    agentsCreated.push_back(agentC);

                }
            }
        }

    private:
        const double maxStrain;
        const double ecmProductionProbability; //0.04; //new ECM agent each 25 hours avg.
        Centerline<3> * centerline;
    };


} // end namespace absmc

#endif
