#ifndef ABSMC_FORCES_3D_H
#define ABSMC_FORCES_3D_H

#include "core/force.h"
#include "core/agentBase.h"
#include "core/centerline.h"
#include "core/trajectorySet.h"
#include "core/mathUtil.h"
#include <util/random.h>
#include <util/agentTypeId.h>
#include <iostream>
#include <limits>
#include <memory>


namespace absmc {

class ZeroUnaryForce3D : public UnaryForce< AgentBase<3> >  {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    point_t calculateForce(agent_t & agent) {
        return point_t(0.,0.,0.);
    }
};


class ZeroBinaryForce3D : public BinaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    ZeroBinaryForce3D() {}

    //calculate interaction force directly
    double calculateForce (const double& firstR, const double& secondR, const double& distance){
        return 0.0;
    }

    //update interaction force between agent1 and agent2 for agent1
    point_t calculateForce(agent_t & agent1, agent_t & agent2) {
        return point_t(0.,0.,0.);
    }
};


/// Composite unary force
/// Forces are calculated in the same order as they are added to the composite rule object.
class CompositeUnaryForce : public UnaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef UnaryForce< AgentBase<3> > force_t;

    void add(force_t* force) {
        forceVector.push_back(force);
    }

    void clear() {
        forceVector.clear();
    }

    virtual void accumulate(agent_t & agent) {
        for (size_t i=0; i<forceVector.size(); i++) {
            (forceVector[i])->accumulate(agent);
        }
    }

    point_t calculateForce(agent_t & agent) {
        point_t force_sum (0.,0.,0.);
        for (size_t i=0; i<forceVector.size(); i++) {
            force_sum += (forceVector[i])->calculateForce(agent);
        }
        return force_sum;
    }

private:
    std::vector<force_t*> forceVector; // vector of different forces included in the composite force
};

/// Type-specific unary force
/// Forces are defined for each agent type
class TypeSpecificUnaryForce : public UnaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef UnaryForce< AgentBase<3> > force_t;

    TypeSpecificUnaryForce(force_t* defaultForce) {
        forceVector.resize(AgentTypeId_MAX, defaultForce);
    }

    void add(force_t* force, AgentTypeId typeId) {
        forceVector[typeId] = force;
    }

    void setDefaultType(AgentTypeId typeId) {
        defaultType = typeId;
    }

    point_t calculateForce(agent_t & agent) {
        return (forceVector[agent.getTypeId()])->calculateForce(agent);
    }

private:
    std::vector<force_t*> forceVector; // vector of different forces; its size is AgentTypeId_MAX;
    AgentTypeId defaultType = tAny;
};

//Gaussian random force; gaussian sigma is equal to input magnitude
class RandomUnaryForce3D : public UnaryForce< AgentBase<3> >  {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;


    RandomUnaryForce3D(double magnitude_, int solverSteps_) :
    magnitude (magnitude_), solverSteps(solverSteps_) {
        randomGen.init();
        curStep = 0;
    }

    point_t calculateForce(agent_t & agent) {
        if (curStep == 0) {
            if (agent.id >= agentForce.size()) {
                agentForce.resize((agent.id+1) * 2);
            }
            //getNextGauss is inefficient, should be rewritten
            double force_x = randomGen.getNextGauss(0, magnitude);
            double force_y = randomGen.getNextGauss(0, magnitude);
            double force_z = randomGen.getNextGauss(0, magnitude);
            agentForce[agent.id] = point_t(force_x, force_y, force_z);
        }

        curStep++;
        if(curStep == solverSteps)
            curStep = 0;

        return agentForce[agent.id];
    }

private:
    double magnitude;
    int solverSteps, curStep;
    util::Random randomGen;
    std::vector<point_t> agentForce;
};


//This force moves synthetic SMCs to the areas with higher concentration of growth factors
//The current realization is very simple: it assumes that:
//1) the concentration gradient is pointed in the direction of the vessel's surface;
//2) the vessel center is at X=0 and the closest surface point is between it and the cell
//a better model would include GF production on the cell-blood interface, its diffusion from there
//and movement in the direction of GF gradient
class GrowthFactorForce3D : public UnaryForce< AgentBase<3> >  {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;


    GrowthFactorForce3D(double magnitude_) :
    magnitude (magnitude_) {}

    point_t calculateForce(agent_t & agent) {
        if((agent.typeId == tSMC3D) &&
           ((dynamic_cast<SMC3D*>(&agent)->getState() == SMC3D::G1) ||
            (dynamic_cast<SMC3D*>(&agent)->getState() == SMC3D::SG2M))) {
            point_t const& pos = agent.pos;
            const point_t center(pos[0], 0, 0); //a slightly better option is to replace this with the centerline of the vessel
            const double distance = metrics.dist(center, pos);
            const double xproject = (center[0]-pos[0]) / distance; //vector from pos to center
            const double yproject = (center[1]-pos[1]) / distance;
            const double zproject = (center[2]-pos[2]) / distance;

            return point_t(xproject * magnitude,
                                   yproject * magnitude,
                                   zproject * magnitude);
        }
        return point_t(0.,0.,0.);
    }

private:
    double magnitude;
    metrics_t metrics;
};


//This force is proportionate to the gradient of WSS, scaling parameter coeff, agents within cutoff used for estimation
//Negative coeff corresponds to migration from high WSS regions to low WSS
class WssGradientForce3D : public UnaryForce< AgentBase<3> >  {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    WssGradientForce3D(double coeff_, double cutoff_) :
    coeff (coeff_), cutoff(cutoff_) {}

    point_t calculateForce(agent_t & agent) {
        if (agent.wssMax == 0.0) {
            return point_t(); //this force should only act on the surface of the vessel
        }

        const size_t nNb = agent.neighbours.size();
        size_t nActive = nNb;
        point_t gradSum(0.0, 0.0, 0.0); //sum of the gradient approximations using k-th neighbour
        double wss = agent.wssMax;
        double cutoffSqr = cutoff * cutoff;
        for (size_t iNb=0; iNb<nNb; iNb++) {
            double wssNb = agent.neighbours[iNb]->wssMax;
            point_t posNb = agent.neighbours[iNb]->pos;
            double distanceSqr = metrics.distSqr(posNb, agent.pos);
            if((wssNb == 0) && (distanceSqr > cutoffSqr)) { //disregard the points which aren't on the surface
                nActive--;
                continue;
            }
            gradSum += (posNb - agent.pos)
                    *((wssNb - wss)/distanceSqr);//finite-diff grad approximation times a unit vector in point-nb direction
        }
        if (nActive == 0)
            return point_t();
        //calculate the average gradient
        gradSum *= (1.0 / (double)nActive);
        //scale
        return gradSum * coeff;
    }

private:
    double coeff, cutoff;
    metrics_t metrics;
};


/// This is applied during force integration to the otherwise immobile stent agents
class UnaryStentDisplacementForce : public UnaryForce< AgentBase<3> >  {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;


    UnaryStentDisplacementForce(double magnitude_, double cY_, double cZ_) :
    magnitude (magnitude_) {
        centerline = new Centerline<3>();
        centerline->resetToAxis(-50., 50., cY_, cZ_, 100); // all reasonable vessels should be inside these coordinates
    }

    UnaryStentDisplacementForce(double magnitude_, Centerline<3> & centerline_) :
    magnitude (magnitude_), centerline(&centerline_) {
    }

    //move the agent away from the center
    point_t calculateForce(agent_t & agent) {
        if ((agent.typeId == tObstacle3D)) {
            point_t pos = agent.getPos();
            // find the closest point on the vessel centerline and move away from it
            const point_t center  = centerline->getClosest(pos);

            // distance from symmetry axis
            const double r = metrics.dist(center, pos);

            // obstacles very close to the central axis are possibly not handled correctly;
            // however these should not occur.
            // dy/y = dr/r
            return r > 0 ? magnitude * (pos-center) / r : point_t(0.,0.,0.);
        }
        return point_t(0.,0.,0.);
    }

private:
    double magnitude;
    Centerline<3> * centerline;
    metrics_t metrics;
};


/// Individual force, kept in a per-agent vector
class IndividualForce3D : public UnaryForce< AgentBase<3> >  {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;


    IndividualForce3D (std::vector<point_t> const& agentForces_) :
    agentForces (agentForces_) {
    }

    IndividualForce3D() {
        agentForces.resize(32000);
    }

    void setForces(std::vector<point_t> const& agentForces_) {
        agentForces = agentForces_;
    }

    void setForce(agent_t const& agent, point_t force) {
        if (agent.getId() >= agentForces.size())
            agentForces.resize(2 * (agent.getId() + 1));
        agentForces[agent.getId()] = force;
    }

    point_t calculateForce(agent_t & agent) {
        return agentForces[agent.getId()];
    }

private:
    std::vector<point_t> agentForces;
};


/// This moves agents along predetermined trajectories
class UnaryTrajectoryForce : public IndividualForce3D {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;


    UnaryTrajectoryForce(double magnitude_, TrajectorySet<3> const& trajSet_, size_t k_ = 16) :
    magnitude (magnitude_), trajSet(trajSet_), k(k_) {
    }

    void updateForce(agent_t * agent, size_t nextStep){
        point_t direction = trajSet.getTrajByInitPosition(agent->getPos0(), nextStep, k);
        setForce (*agent, direction * magnitude);
    }


private:
    double magnitude;
    TrajectorySet<3> trajSet;
    size_t k; //number of neighbours to sample
};


/// This is applied for inflation testing
/// Relies on precalculated surface normals
class UnaryInflationForce : public UnaryForce< AgentBase<3> >  {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;


    UnaryInflationForce(double magnitude_, double cY_, double cZ_) :
    magnitude (magnitude_) {
        centerline = new Centerline<3>();
        centerline->resetToAxis(-50., 50., cY_, cZ_, 100); // all reasonable vessels should be inside these coordinates
    }

    UnaryInflationForce(double magnitude_, Centerline<3> & centerline_) :
    magnitude (magnitude_), centerline(&centerline_) {
    }

    void setMagnitude(double magnitude_) {
        magnitude = magnitude_;
    }

    //move the agent away from the center by pushing it along the surface normal
    point_t calculateForce(agent_t & agent) {
        point_t pos = agent.getPos();
        // apply to all cells on the inner surface
        if((agent.isSurface()) && centerline->isFacingCenterline(pos, *agent.getSurfaceNormal())) {
            // Apply pressure in the direction opposite to the surface normal
            return (*agent.getSurfaceNormal()) * (- magnitude);
        }
        return point_t(0.,0.,0.);
    }

private:
    double magnitude;
    Centerline<3> * centerline;
};


/// generic pairwise spherically symmetric 3D interaction force
class BinaryForce3D : public BinaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    //Force > 0 means attraction
    virtual double calculateForce(const double& firstR, const double& secondR, const double& distance) = 0;

    virtual point_t calculateForce(agent_t & agent1, agent_t & agent2) {

        // distance (real dimensions)
        const double distance = metrics.dist(agent1.pos, agent2.pos);
        const double force = calculateForce(agent1.r, agent2.r, distance);

        if(force == 0.)
            return point_t();

        const double xproject = (agent2.pos[0]-agent1.pos[0]) / distance;
        const double yproject = (agent2.pos[1]-agent1.pos[1]) / distance;
        const double zproject = (agent2.pos[2]-agent1.pos[2]) / distance;

        return point_t(xproject * force,
                        yproject * force,
                        zproject * force);
    }

protected:
    metrics_t metrics;
};


/// Composite binary force
/// Forces are calculated in the same order as they are added to the composite rule object.
class CompositeBinaryForce : public BinaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef BinaryForce< AgentBase<3> > force_t;

    void add(force_t* force) {
        forceVector.push_back(force);
    }

    void clear() {
        forceVector.clear();
    }

    point_t calculateForce(agent_t & agent1, agent_t & agent2) {
        point_t force_sum (0.,0.,0.);
        for (size_t i=0; i<forceVector.size(); i++) {
            force_sum += (forceVector[i])->calculateForce(agent1, agent2);
        }
        return force_sum;
    }

private:
    std::vector<force_t*> forceVector; // vector of different forces included in the composite force
};


/// Type-specific binary force
/// Forces are defined for symmetric pairs of agent types
class TypeSpecificBinaryForce : public BinaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef BinaryForce< AgentBase<3> > force_t;

    TypeSpecificBinaryForce(force_t* defaultForce) {
        forceVector.resize(AgentTypeId_MAX * AgentTypeId_MAX, defaultForce);
    }

    // forces are added symmetrically for both (type1,type2) and (type2,type1)
    void add(force_t* force, AgentTypeId typeId1, AgentTypeId typeId2) {
        forceVector[typeId1 * AgentTypeId_MAX + typeId2] =
        forceVector[typeId2 * AgentTypeId_MAX + typeId1] = force;
    }

    void setDefaultType(AgentTypeId typeId) {
        defaultType = typeId;
    }

    point_t calculateForce(agent_t & agent1, agent_t & agent2) {
        return (forceVector[agent1.getTypeId() * AgentTypeId_MAX + agent2.getTypeId()])->calculateForce(agent1, agent2);
    }

    /// Returns interaction force for the default cell type
    double calculateForce(const double& firstR, const double& secondR, const double& distance) {
        return (forceVector[defaultType * AgentTypeId_MAX + defaultType])->calculateForce(firstR, secondR, distance);
    }

private:
    std::vector<force_t*> forceVector; // vector of different pairwise forces; its size is AgentTypeId_MAX * AgentTypeId_MAX;
    AgentTypeId defaultType = tAny;
};


class NeohookeanRepulsionForce3D : public BinaryForce3D {
public:

    /// Constructor; characteristic length scale L in real units.
    NeohookeanRepulsionForce3D(double L_) :
    L(L_) { }

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE
    //Force > 0 means attraction
    double calculateForce (const double& firstR, const double& secondR, const double& distance){

        //INDENTATION OF SOFT MATTER BEYOND THE HERTZIAN REGIME
        //NEOHOOKEAN
        const double Rsum = firstR + secondR;
        if (distance > Rsum){
            return 0;
        }
        else {
            const double identdepth = fabs(Rsum - distance);
            const double Rad = firstR*secondR/Rsum;//1.0/(1.0/R1 + 1.0/R2); //average Radius of the sphere
            const double Rsqr = Rad*Rad;
            //const double Rcub = Rsqr*Rad;

            const double asqr = Rad*identdepth;
            const double a = sqrt(asqr); //contact radius, taken from identation depth via sqrt(R*identdepth) or identdepth^0.498 * R^(0.503-3.97e-6*indetdepth),see ref Lin et al Spherical identation of soft matter
            // we could also try to measure contact radius directly as it is done in isr2d

            const double acub = a*asqr;
            const double afou = asqr*asqr;
            const double afif = acub*asqr;
            //const double firstpart = (afif - 15.0*Rad*afou + 75.0*Rsqr*acub) / (5.0*Rad*asqr - 50.0*Rsqr*a + 125.0*Rcub) ;
            const double firstpart = (16.0*afif - 36.0*math::pi*Rad*afou + 27.0*math::pi*math::pi*Rsqr*acub) / (4.0*a-3.0*math::pi*Rad)/(4.0*a-3.0*math::pi*Rad) ; // mm^3
            //double secondpart = (acub - 15.0*Rad*asqr) / (25.0*Rsqr*a - 125.0*Rcub);
            const double B = 0.2; //MPa  //B in terms of Young
            const double f12rad = - 8.0*B * firstpart/ 3.0 / Rad /** exp(bsmall*secondpart)*/ ; //MPa * mm^2 = N
            return f12rad;
        }
    }

private:
    double L;
};

/// Sigmoid SMC attraction force; attempt to replicate macroscopic behavoiur without bond breaking
class SMCAttractionForce3D : public BinaryForce3D {
public:

    /// Constructor; characteristic length scale L in real units.
    SMCAttractionForce3D (double L_) :
    L(L_) { }

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE
    //Force > 0 means attraction
    double calculateForce (const double& firstR, const double& secondR, const double& distance){

        //INDENTATION OF SOFT MATTER BEYOND THE HERTZIAN REGIME
        //NEOHOOKEAN
        const double Rsum = firstR + secondR;

        if (distance <= Rsum) {
            return 0;
        }
        else {

            //linear interaction force
            /*const double linconst = 0.004; // N/mm; based on data for media in Holzapfel et al., 2005 // alternative: 0.06
            const double maxdist = 1.3;
            const double identdepth = fabs(Rsum - distance);
            if (distance > maxdist * Rsum){
                f12rad =  linconst * (maxdist - 1) * Rsum;
            }
            else {
                f12rad =  linconst * identdepth;
            }*/

            // sigmoidal interaction force
            const double strain = distance / Rsum;
            const double c1 = 0.001 * 150;
            const double c2 = 15.;
            const double c3 = 11.5;
            const double c4 = c1*(c2-c3) / sqrt((c2-c3)*(c2-c3) + 1); //shift value at strain = 1 to zero
            const double t  = c3 * strain - c2;
            const double sigmoid = c1 * t / sqrt(1 + t * t) + c4; // interaction force based on circumferential stress-strain curve for media tissue from Holzapfel et al., 2005; force in MPa = N/mm^2
            // real tissue breaks at ~1.4 strain, this function starts to flatten out at this point
            // adjusted for the number of contacts in a hexagonal lattice per mm^2
            const double f12rad = sigmoid * (Rsum * Rsum) * 0.3536; // 1 / (2 * sqrt(2)) = 0.3536  // force in N
            return f12rad;
        }
    }

private:
    double L;
};


/// Linear force calibrated for the average stress-atrain relationship used for animal models by POLIMI
class LinearAttractionForce3D : public BinaryForce3D {
public:

    /// Constructor;
    LinearAttractionForce3D (double c1) :
    c1(c1) { }

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE
    //Force > 0 means attraction
    double calculateForce (const double& firstR, const double& secondR, const double& distance){

        //INDENTATION OF SOFT MATTER BEYOND THE HERTZIAN REGIME
        //NEOHOOKEAN
        const double Rsum = firstR + secondR;

        if (distance <= Rsum) {
            return 0;
        }
        else {

            //linear interaction force
            /*const double linconst = 0.004; // N/mm; based on data for media in Holzapfel et al., 2005 // alternative: 0.06
            const double maxdist = 1.3;
            if (distance > maxdist * Rsum){
                f12rad =  linconst * (maxdist - 1) * Rsum;
            }
            else {
                f12rad =  linconst * identdepth;
            }*/

            const double identdepth = fabs(Rsum - distance);
            const double macro = c1 * identdepth / Rsum;
            // adjusted for the number of contacts in a hexagonal lattice per mm^2
            return macro * (Rsum * Rsum) * 0.3536; // 1 / (2 * sqrt(2)) = 0.3536
        }
    }
private:
    double c1;
};


// f = c1 * l^4 + c2 * l
class Linear4thPowerAttractionForce3D : public BinaryForce3D {
public:

    /// Constructor;
    Linear4thPowerAttractionForce3D (double L_) :
    L(L_) { }

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE
    //Force > 0 means attraction
    double calculateForce (const double& firstR, const double& secondR, const double& distance){

        const double Rsum = firstR + secondR;

        if (distance <= Rsum) {
            return 0;
        }
        else {

            //linear interaction force
            /*const double linconst = 0.004; // N/mm; based on data for media in Holzapfel et al., 2005 // alternative: 0.06
            const double maxdist = 1.3;
            if (distance > maxdist * Rsum){
                f12rad =  linconst * (maxdist - 1) * Rsum;
            }
            else {
                f12rad =  linconst * identdepth;
            }*/

            const double identdepth = fabs(Rsum - distance);
            const double c1 = 0.001 * 12000;
            const double c2 = 0.001 * 50; // linear component for small deformations
            const double macro = c1 * (identdepth / Rsum) * (identdepth / Rsum) * (identdepth / Rsum) * (identdepth / Rsum) + c2 * (identdepth / Rsum);
            // adjusted for the number of contacts in a hexagonal lattice per mm^2
            return macro * (Rsum * Rsum) * 0.3536; // 1 / (2 * sqrt(2)) = 0.3536
        }
    }
private:
    double L;
};


/// 6th power polynomial. Used for media and adventitia in POLIMI three-layer model
class Poly6thPowerAttractionForce3D : public BinaryForce3D {
public:

    /// Constructor; coefficients from 1th to 6th power
    /// the scaling coefficient is empirical, introduced to get the desired macroscopic behaviour
    Poly6thPowerAttractionForce3D (double c1, double c2, double c3, double c4, double c5, double c6) :
        c1(c1), c2(c2), c3(c3), c4(c4), c5(c5), c6(c6)
    { }

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE [mm, N, MPa]
    //Force > 0 means attraction
    double calculateForce (const double& firstR, const double& secondR, const double& distance){

        const double Rsum = firstR + secondR;

        if (distance <= Rsum) {
            return 0;
        }
        else {
            const double identdepth = fabs(Rsum - distance);

            const double strain = identdepth / Rsum;
            const double strain2 = strain * strain;
            const double strain3 = strain * strain2;
            const double strain4 = strain * strain3;
            const double strain5 = strain * strain4;
            const double strain6 = strain * strain5;

            const double macro = c6 * strain6 + c5 * strain5 + c4 * strain4 + c3 * strain3 + c2 * strain2 + c1 * strain;
            // scaling coefficient from macro to micro force
            // adjusted for the number of contacts in a hexagonal lattice per mm^2
            // 1 / (2 * sqrt(2)) = 0.3536

            return macro * (Rsum * Rsum);
        }
    }

private:
    double c1, c2, c3, c4, c5, c6;
};


/// BELOW THIS POINT BE DRAGONS
/// Old forces, not updated anymore

/// Bilinear plus Morse force with unlimited interaction range.
/// Deprecated; if you see it used, split it into bond + repulsion forces
class Bilinear3D : public BinaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    /// Constructor; characteristic length scale L in real units.
    Bilinear3D(double L_) :
    L(L_) { }

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE
    //Force > 0 means attraction
    double calculateForce (const double& firstR, const double& secondR, const double& distance){

                // ALTERNATIVE 1
                //        train1 = fabs(train1);
                //        f12 = 0.25*118.34*exp(13.2*train1)*0.0942477796*aux/25.66;        //dividim per 4 a veure si va millor

                // ALTERNATIVE 2
                    /*

                    if (0.021<r) {f12 = 17.13*(r-0.03);} else {
                        if (r>=0.039) {f12 = 97.66*(r-(0.039-17.13*0.009/97.66)); } else {
                            f12 = 97.66*(r-(0.021-17.13*0.009/97.66));
                        }
                    }


                    if (0.0279<r && r<0.039) {f12 = 1.713*(r-0.03);} else {
                        if (r<0.0279 ) {f12 = -U1OverR1 * exp(-(r+0.002)*OneOverLR1) + U2OverR2 * exp(-(r+0.002)*OneOverLR2);} else {
                            if (r>=0.039 && r<0.0476) {f12 = 36.0894*(r-(0.039-1.713*0.009/36.0894)); } else {
                                f12 = -U1OverR1 * exp(-r*OneOverLR1) + U2OverR2 * exp(-r*OneOverLR2);
                            }
                        }
                    }
                    */


                ////        f12 = (3.14159*0.015/12.0)*(xdist-xdist0)/L;//118.34*exp(13.2*train)*0.0942477796*(r/L-2.0)/25.66;

        // ALTERNATIVE 3, INDENTATION OF SOFT MATTER BEYOND THE HERTZIAN REGIME

        //NEOHOOKEAN
        const double Rsum = firstR + secondR;
        const double identdepth = fabs(Rsum - distance);

        double f12rad;

        if (distance > Rsum){

            //linear interaction force
            /*const double linconst = 0.004; // N/mm; based on data for media in Holzapfel et al., 2005 // alternative: 0.06
            const double maxdist = 1.3;
            if (distance > maxdist * Rsum){
                f12rad =  linconst * (maxdist - 1) * Rsum;
            }
            else {
                f12rad =  linconst * identdepth;
            }*/

            // sigmoidal interaction force
            const double strain = distance / Rsum;
            const double c1 = 0.001 * 150;
            const double c2 = 15.;
            const double c3 = 11.5;
            const double c4 = c1*(c2-c3) / sqrt((c2-c3)*(c2-c3) + 1); //shift value at strain = 1 to zero
            const double t  = c3 * strain - c2;
            const double sigmoid = c1 * t / sqrt(1 + t * t) + c4; // interaction force based on circumferential stress-strain curve for media tissue from Holzapfel et al., 2005; force in MPa = N/mm^2
            // real tissue breaks at ~1.4 strain, this function starts to flatten out at this point
            // adjusted for the number of contacts in a hexagonal lattice per mm^2
            f12rad = sigmoid * (Rsum * Rsum) * 0.3536; // 1 / (2 * sqrt(2)) = 0.3536  // force in N

        } else {
            const double Rad = firstR*secondR/Rsum;//1.0/(1.0/R1 + 1.0/R2); //average Radius of the sphere
            const double Rsqr = Rad*Rad;
            //const double Rcub = Rsqr*Rad;

            const double asqr = Rad*identdepth;
            const double a = sqrt(asqr); //contact radius, taken from identation depth via sqrt(R*identdepth) or identdepth^0.498 * R^(0.503-3.97e-6*indetdepth),see ref Lin et al Spherical identation of soft matter

            // we could also try to measure contact radius directly as it is done in isr2d
            const double acub = a*asqr;
            const double afou = asqr*asqr;
            const double afif = acub*asqr;
            //const double firstpart = (afif - 15.0*Rad*afou + 75.0*Rsqr*acub) / (5.0*Rad*asqr - 50.0*Rsqr*a + 125.0*Rcub) ;
            const double firstpart = (16.0*afif - 36.0*math::pi*Rad*afou + 27.0*math::pi*math::pi*Rsqr*acub) / (4.0*a-3.0*math::pi*Rad)/(4.0*a-3.0*math::pi*Rad) ; // mm^3

    //        double secondpart = (acub - 15.0*Rad*asqr) / (25.0*Rsqr*a - 125.0*Rcub);
            const double B = 0.2; //MPa //0.0282385035;//1.0/8.3439;//0.5*19.59; //kPa  //B in terms of young, looking for efficiency
            f12rad = - 8.0*B * firstpart/ 3.0 / Rad /** exp(bsmall*secondpart)*/ ; //MPa * mm^2 = N
        //    f12rad =  -4.0*acub*youngeq/(0.015*3.0*Rad);
        }

        //        if (tEEL3D == agent2.getTypeId()) {
        //            f12rad = 0.5*f12rad;
        //        }
        //        agent1.setStrain(lambda);

        return f12rad;
    }

    point_t calculateForce(agent_t & agent1, agent_t & agent2) {
        point_t const& pos1 = agent1.pos;
        point_t const& pos2 = agent2.pos;

        // distance (real dimensions)
        const double distance = metrics.dist(pos1, pos2);
        const double R1 = agent1.r;
        const double R2 = agent2.r;

        const double force = calculateForce(R1, R2, distance);
        //extra attraction if membranes are involved
        // if ((agent1.typeId == tIEL3D) && (agent2.typeId == tIEL3D) && (force > 0)) {
        //     force = force * 2.0;
        // }
        //
        // if ((agent1.typeId == tEEL3D) && (agent2.typeId == tEEL3D) && (force > 0)) {
        //     force = force * 2.0;
        // }

        const double xproject = (pos2[0]-pos1[0]) / distance;
        const double yproject = (pos2[1]-pos1[1]) / distance;
        const double zproject = (pos2[2]-pos1[2]) / distance;

        return point_t(xproject * force,
                        yproject * force,
                        zproject * force);
    }

    void accumulate(agent_t & agent1, agent_t & agent2) {

        point_t const& pos1 = agent1.pos;
        point_t const& pos2 = agent2.pos;

        // distance (real dimensions)
        const double distance = metrics.dist(pos1, pos2);
        const double R1 = agent1.r;
        const double R2 = agent2.r;

        const double force = calculateForce(R1, R2, distance);

        //extra attraction if membranes are involved
        // if ((agent1.typeId == tIEL3D) && (agent2.typeId == tIEL3D) && (force > 0)) {
        //     force = force * 2.0;
        // }
        //
        // if ((agent1.typeId == tEEL3D) && (agent2.typeId == tEEL3D) && (force > 0)) {
        //     force = force * 2.0;
        // }

        const double xproject = (pos2[0]-pos1[0]) / distance;
        const double yproject = (pos2[1]-pos1[1]) / distance;
        const double zproject = (pos2[2]-pos1[2]) / distance;

        agent1.force += point_t(xproject * force, //2.0 for force anisotropy (probably not a very good idea)
                                yproject * force,
                                zproject * force);
/*
        const double l1 = R1 + 0.5*(distance - R1 - R2);
        const double f12radeq = l1*force; //virial

        agent1.addStress(point_t(xproject * f12radeq,   //not needed on intermediate steps
                                 yproject * f12radeq,       //used in the necrosis rule
                                 zproject * f12radeq) , R1*R1*R1);
        const double f_times_l_over_v = 0.75 * f12radeq/(3.1415926 * R1*R1*R1); // minus to make stress matrix components (and pressure) positive when forces are repulsive
        agent1.addStressMatrix(f_times_l_over_v * xproject * xproject, // virial stress matrix, eq. to macroscopic Cauchy stress
                                f_times_l_over_v * yproject * xproject,                    // elements are f_x * delta-r_x etc.
                                f_times_l_over_v * zproject * xproject,
                                f_times_l_over_v * yproject * yproject,
                                f_times_l_over_v * zproject * zproject,
                                f_times_l_over_v * zproject * yproject); */ // 00, 10, 20, 11, 22, 21
    }

private:
    const double L;
    metrics_t metrics;
};


/// Neohookean repulsion + Lennard-Jones attraction with cutoff radius
class NeohookeanLJ3D : public BinaryForce< AgentBase<3> > {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    /// Constructor; characteristic length scale L_ and cutoff distance cutoff_ in real units.
    NeohookeanLJ3D (double L_, double cutoff_) :
    L(L_), cutoff(cutoff_), epsilon(0.0000002) { } //06 //epsilon(0.0000001)

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE (distance in mm)
    //Force > 0 means attraction
    double calculateForce (double firstR, double secondR, double distance){
        //NEOHOOKEAN
        const double Rsum = firstR + secondR;

        const double identdepth = fabs(Rsum - distance);

        double f12rad;

        if (distance > Rsum){
            if (distance > cutoff) {
                f12rad = 0;
            }
            //LJ attraction
            else {
                const double sigmaByR = Rsum / distance;
                const double epsilonByR = epsilon / distance;
                const double sigmaByRcub = sigmaByR * sigmaByR * sigmaByR;
                const double sigmaByRsix = sigmaByRcub * sigmaByRcub;
                const double sigmaByR12th = sigmaByRsix * sigmaByRsix;
                f12rad = - 48 * epsilonByR * (sigmaByR12th - sigmaByRsix);
                //f12rad =  0.01*identdepth;
            }
        } else {
            const double Rad = firstR*secondR/Rsum;//1.0/(1.0/R1 + 1.0/R2); //average Radius of the sphere
            const double Rsqr = Rad*Rad;
            //const double Rcub = Rsqr*Rad;

            const double asqr = Rad*identdepth;
            const double a = sqrt(asqr); //contact radius, taken from identation depth via sqrt(R*identdepth) or identdepth^0.498 * R^(0.503-3.97e-6*indetdepth),see ref Lin et al Spherical identation of soft matter

            // we could also try to measure contact radius directly as it is done in isr2d
            const double acub = a*asqr;
            const double afou = asqr*asqr;
            const double afif = acub*asqr;
            const double firstpart = (16.0*afif - 36.0*math::pi*Rad*afou + 27.0*math::pi*math::pi*Rsqr*acub) / (4.0*a-3.0*math::pi*Rad)/(4.0*a-3.0*math::pi*Rad) ;

            const double B = 0.1; //MPa //0.0282385035;//1.0/8.3439;//0.5*19.59; //kPa  //B in terms of young, looking for efficiency
            f12rad = - 8.0*B * firstpart/ 3.0 / Rad /** exp(bsmall*secondpart)*/ ;
        }
        //if (std::isnan(f12rad)) f12rad = 0.0;
        return f12rad;
    }

    void accumulate(agent_t & agent1, agent_t & agent2) {

        //case 1: both agents are mobile and the forces have already been calculated
        if (agent2.isMobile && agent1.id > agent2.id)
            return;

        //cases 2 and 3: forces have not been calculated
        point_t const& pos1 = agent1.pos;
        point_t const& pos2 = agent2.pos;

        // distance (real dimensions)
        const double distance = metrics.dist(pos1, pos2);
        const double R1 = agent1.r;
        const double R2 = agent2.r;

        //DEBUG
        /*if ((agent1.typeId == tObstacle3D) || (agent2.typeId == tObstacle3D)) {
            distance = pos2[1]-pos1[1];
        }*/

        double force = calculateForce(R1, R2, distance);

        //extra attraction if membranes are involved
        if ((agent1.typeId == tIEL3D) && (agent2.typeId == tIEL3D) && (force > 0)) {
            force = force * 2.0;
        }

        if ((agent1.typeId == tEEL3D) && (agent2.typeId == tEEL3D) && (force > 0)) {
            force = force * 2.0;
        }
        //if ((agent1.typeId == tEndo3D) && (agent2.typeId == tEndo3D) && (force > 0)) {
        //    force = 0.0;
        //}

        //walls aren't squishy
        //if (((agent1.typeId == tObstacle3D) || (agent2.typeId == tObstacle3D)) && (force < 0)) {
        //    force = force * 5.0;
        //}

        //Also, increased adhesion to walls
        //if (((agent1.typeId == tObstacle3D) || (agent2.typeId == tObstacle3D)) && (force > 0)) {
        //    force = force * 2.0;
        //}


        const double xproject = (pos2[0]-pos1[0]) / distance;
        const double yproject = (pos2[1]-pos1[1]) / distance;
        const double zproject = (pos2[2]-pos1[2]) / distance;

        //DEBUG
        /*if ((agent1.typeId == tObstacle3D) || (agent2.typeId == tObstacle3D)) {
            xproject = 0;
            yproject = 1;
            zproject = 0;
        }*/

        const double l1 = R1 + 0.5*(distance - R1 - R2);
        const double f12radeq = l1*force; //virial

        //case 2: if agent 2 is immobile, only add the force to agent1
        agent1.force += point_t(xproject * force, //2.0 for force anisotropy (probably not a very good idea)
                                yproject * force,
                                zproject * force);
        agent1.addStress(point_t(xproject * f12radeq,    //not needed on intermediate steps
                                 yproject * f12radeq,
                                 zproject * f12radeq) , R1*R1*R1);
        //case 3: if agent 2 is mobile, and the force has not been added for it, add the opposite force
        if (agent2.isMobile) {
            const double l2 = R2 + 0.5*(distance - R1 - R2);
            const double f21radeq = l2*force; //virial
            agent2.force += point_t(- xproject * force, //2.0 for force anisotropy (probably not a very good idea)
                                    - yproject * force,
                                    - zproject * force);
            agent2.addStress(point_t(xproject * f21radeq,    //not needed on intermediate steps
                                     yproject * f21radeq,
                                     zproject * f21radeq) , R1*R1*R1);
        }
    }

//    void calculatedensity(agent_t & agent1, agent_t const& agent2) {

//        point_t const& pos1 = agent1.getPos();
//        point_t const& pos2 = agent2.getPos();

//        const double xdist = fabs(pos2[0]-pos1[0]); //inter-particle distance

//        double dist = 4.0*L;
//        double z = xdist/dist;
//        double W = 0.0;

//        // JB NOTE: z==1 and z >= 2 are not handled, is this on purpose?
//        if (z<1) {
//            W = 2.0*(1.0-1.5*z*z+0.75*z*z*z)/(3.0*dist);
//        }
//        else if (1<z && z<2) {
//            W = (2.0-z)*(2.0-z)*(2.0-z)/(6.0*dist);
//        }

//        agent1.addDensity(W);
//    }

private:
    const double L;
    metrics_t metrics;
    double cutoff;
    double epsilon;
};


/// worm-like chain potential + hertzian repulsion. Ref: N. Melnikova et al., 2017
class Wmlc3D : public BinaryForce< AgentBase<3> > {
private:
    double const lmax_c, hertz_x, lmax_x, hertz_c,
        wmlc_x, wmlc_c;

public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    /// Constructor; characteristic length scale L in real units.
    Wmlc3D(double L_) :
        lmax_c(1.27), hertz_x(0.43*0.02), lmax_x(1.10),
        hertz_c ((1.4+3*0.1)*1.52*16*0.02),
        wmlc_x (0.00000007*(0.7+0.1*1)*1.46*0.02),
        wmlc_c ((0.00000025+0.00000002*3)*0.9*16*0.02),
        L(L_)
    {

    }

    //calculate interaction force between two agents
    //EVERYTHING IN REAL UNITS HERE
    double calculateForce (const double& firstR, const double& secondR, const double& distance, const double& cosphi, const bool& useWmlc = true){
        //NEOHOOKEAN
        const double Rsum = firstR + secondR;
        const double identdepth = fabs(Rsum - distance);
        double f12rad;

         //wmlc+hertz model
        const double lmax = (firstR+secondR)*(lmax_x+(lmax_c-lmax_x)*cosphi);
        const double Rad = firstR*secondR/Rsum;//1.0/(1.0/R1 + 1.0/R2); //average Radius of the sphere
        const double zz=distance/lmax;
        const double kbT=wmlc_x+(wmlc_c-wmlc_x)*cosphi;
        const double hertz_coeff = hertz_x - (hertz_x-hertz_c)*cosphi;
        const double popravka=0.036/Rad;

        if (distance <= Rsum){
            const double Rsqr = Rad*Rad;
            const double asqr = Rad*identdepth;
            const double a = sqrt(asqr); //contact radius, taken from identation depth via sqrt(R*identdepth) or identdepth^0.498 * R^(0.503-3.97e-6*indetdepth),see ref Lin et al Spherical identation of soft matter

            // we could also try to measure contact radius directly as it is done in isr2d
            const double acub = a*asqr;
            const double afou = asqr*asqr;
            const double afif = acub*asqr;
            const double znam = (4.0*a-3.0*math::pi*Rad)*(4.0*a-3.0*math::pi*Rad);
            const double firstpart = (16.0*afif - 36.0*math::pi*Rad*afou + 27.0*math::pi*math::pi*Rsqr*acub) / znam;
            const double B = 0.0282385035; //MPa
            //const double B = 0.1; //MPa //0.0282385035;//1.0/8.3439;//0.5*19.59; //kPa  //B in terms of young, looking for efficiency
            f12rad = - 8.0 * B * firstpart / 3.0 / Rad;
            f12rad = f12rad * hertz_coeff;
        } else {
            f12rad=0.;
        }

        double fwmlc=(zz<0.99)?(kbT/lmax/(popravka*popravka*popravka)*(6.0*zz+(3.0*zz*zz-2.0*zz*zz*zz)/(1-zz)/(1-zz))):(kbT*10000.0/lmax/(popravka*popravka*popravka));

        return useWmlc ? f12rad+fwmlc : f12rad;
        //end of wmlc+hertz model    */
    }

    //return the interaction force for "average" configuration. Used to determine the characteristic force.
    double calculateForce (const double& firstR, const double& secondR, const double& distance)
    {
        return calculateForce (firstR, secondR, distance, 0.5);
    }

    point_t calculateForce(agent_t & agent1, agent_t & agent2) {

        point_t const& pos1 = agent1.pos;
        point_t const& pos2 = agent2.pos;

        // distance (real dimensions)
        const double distance = metrics.dist(pos1, pos2);
        const double R1 = agent1.r;
        const double R2 = agent2.r;

        //here the initial angle is used
        const double cosfi = calculateAngleCos(agent1.pos0, agent2.pos0); //angle between spring and circumferential direction

        double force;
        if(agent1.getTypeId() == tObstacle3D || agent2.getTypeId() == tObstacle3D)
            force = calculateForce(R1, R2, distance, cosfi, false);
        else
            force = calculateForce(R1, R2, distance, cosfi);

        const double xproject = (pos2[0]-pos1[0]) / distance;
        const double yproject = (pos2[1]-pos1[1]) / distance;
        const double zproject = (pos2[2]-pos1[2]) / distance;

        return point_t(xproject * force,
                                yproject * force,
                                zproject * force);
    }

//    void calculatedensity(agent_t & agent1, agent_t const& agent2) {

//        point_t const& pos1 = agent1.getPos();
//        point_t const& pos2 = agent2.getPos();

//        const double xdist = fabs(pos2[0]-pos1[0]); //inter-particle distance

//        double dist = 4.0*L;
//        double z = xdist/dist;
//        double W = 0.0;

//        // JB NOTE: z==1 and z >= 2 are not handled, is this on purpose?
//        if (z<1) {
//            W = 2.0*(1.0-1.5*z*z+0.75*z*z*z)/(3.0*dist);
//        }
//        else if (1<z && z<2) {
//            W = (2.0-z)*(2.0-z)*(2.0-z)/(6.0*dist);
//        }

//        agent1.addDensity(W);
//    }

private:

    //calculate the angle cosine between R1R2 and circumferential direction. 0 < angle < pi
    inline double calculateAngleCos(point_t const& pos1, point_t const& pos2)
    {
        //assuming pos1 and pos2 are close, we use pos1 to determine the circumferential direction
        //that is, (0, pos1[1], pos1[2]), rotated pi/2 in YZ plane.
        //or, (0, -pos1[2], pos1[1]); then we calculate a normalized vector in that direction
        const double circlength = sqrt((pos1[1]*pos1[1])+(pos1[2]*pos1[2]));
        point_t circdir (0, -pos1[2]/circlength, pos1[1]/circlength);
        //next, the length of the projection of R1R2 on circdir.
        const double circproj = abs( (pos2[1]-pos1[1])*circdir[1] + (pos2[2]-pos1[2])*circdir[2] );
        const double distance = metrics.dist(pos1, pos2);

        //if ((!((circproj/distance)<=1)) && (!((circproj/distance)>=0))) {
        //    std::cout << "OH NO " << circproj/distance << " " << pos1[0] << " " << pos2[0] << endl;
        //}

        return circproj/distance; //angle between spring and circumferential direction
    }

    const double L;
    metrics_t metrics;
};

/// Morse force with unlimited interaction range.
/// NOTE: new version, handles non-unity length scales.
class MorseForce3D {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    /// Constructor; U1,U2,r1,r2,rCutOff to be specified in non-dimensional system, characteristic length scale L in real units.
    MorseForce3D(double U1_, double r1_, double U2_, double r2_, double L_) :
        U1(U1_), r1(r1_), U2(U2_), r2(r2_), L(L_),
        U1OverR1(U1/r1), U2OverR2(U2/r2), OneOverLR1(1.0/(L*r1)), OneOverLR2(1.0/(L*r2))
        { }

    void accumulate(agent_t & agent1, agent_t & agent2) {

        // distance (real dimensions)
        const double r = metrics.dist(agent1.getPos(), agent2.getPos());
        // Morse force (non-dimensional)
        const double f12 = -U1OverR1 * exp(-r*OneOverLR1) + U2OverR2 * exp(-r*OneOverLR2);

        #ifdef ABSMC_EXCESSIVE_CHECKING
        if (isnan(f12)) std::cerr << "morse_force: isnan(f12)\n";
        if (isnan(1.0/r)) std::cerr << "morse_force: isnan(1/r)\n";
        #endif

        point_t const& pos1 = agent1.getPos();
        point_t const& pos2 = agent2.getPos();

        agent1.addForce(point_t((pos2[0]-pos1[0]) / r * f12,
                                (pos2[1]-pos1[1]) / r * f12,
                                (pos2[2]-pos1[2]) / r * f12) );
    }

    void accumulateforce(agent_t & agent1, agent_t const& agent2) { }
    void accumulatestress(agent_t & agent1, agent_t const& agent2) { }
    void calculatedensity(agent_t & agent1, agent_t const& agent2) { }
private:
    const double U1, r1, U2, r2, L;
    const double U1OverR1, U2OverR2, OneOverLR1, OneOverLR2;
    metrics_t metrics;
};


/// Morse force with cut-off range.
/// NOTE: new version, handles non-unity length scales.
class CutOffMorseForce3D {
public:
    typedef AgentBase<3> agent_t;
    typedef Point<agent_t::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    /// Constructor; U1,U2,r1,r2,rCutOff to be specified in non-dimensional system, characteristic length scale L in real units.
    CutOffMorseForce3D(double U1_, double r1_, double U2_, double r2_, double rCutOff_, double L_) :
        U1(U1_), r1(r1_), U2(U2_), r2(r2_), rCutOff(rCutOff_), L(L_),
        U1OverR1(U1/r1), U2OverR2(U2/r2), OneOverLR1(1.0/(L*r1)), OneOverLR2(1.0/(L*r2)), rCutOffTimesL(rCutOff*L)
        { }

    void accumulate(agent_t & agent1, agent_t & agent2) {

        // distance (real dimensions)
        const double r = metrics.dist(agent1.getPos(), agent2.getPos());
        if (r>rCutOffTimesL) return;

        // Morse force (non-dimensional)
        const double f12 = -U1OverR1 * exp(-r*OneOverLR1) + U2OverR2 * exp(-r*OneOverLR2);

        #ifdef ABSMC_EXCESSIVE_CHECKING
        if (isnan(f12)) std::cerr << "morse_force: isnan(f12)\n";
        if (isnan(1.0/r)) std::cerr << "morse_force: isnan(1/r)\n";
        #endif

        point_t const& pos1 = agent1.getPos();
        point_t const& pos2 = agent2.getPos();

        agent1.addForce(point_t((pos2[0]-pos1[0]) / r * f12,
                                (pos2[1]-pos1[1]) / r * f12,
                                (pos2[2]-pos1[2]) / r * f12) );
    }

    void accumulateforce(agent_t & agent1, agent_t const& agent2) { }
    void accumulatestress(agent_t & agent1, agent_t const& agent2) { }
    void calculatedensity(agent_t & agent1, agent_t const& agent2) { }
private:
    const double U1, r1, U2, r2, rCutOff, L;
    const double U1OverR1, U2OverR2, OneOverLR1, OneOverLR2, rCutOffTimesL;
    metrics_t metrics;
};

} // end namespace absmc

#endif
