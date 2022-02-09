#ifndef ABSMC_RUNGE_KUTTA_INTEGRATOR_H
#define ABSMC_RUNGE_KUTTA_INTEGRATOR_H

#include <vector>
#include <iostream>

#include "core/geometry.h"
#include "core/integrator.h"
#include "core/neighbourDetector.h"
#include "core/force.h"
#include <util/simpleStatistics.h>

namespace absmc {

/// Integrator using a RK4 scheme with adaptive step size
/// to integrate the overdamped equation of motion.
template<class Agent>
class RungeKuttaIntegratorMPI : public Integrator<Agent> {
public:
    typedef Point<Agent::nDim, double> point_t;
    typedef BinaryForce < Agent > binForce;
    typedef UnaryForce < Agent > unForce;

    /// Constructor; characteristic length scale L and maxDispl to be specified in real units, maxTimeStep in non-dimensional system.
    RungeKuttaIntegratorMPI(unForce & unaryForces_, binForce & binaryForces_, BondForce <Agent> & bondForces_,
                            double L_,
                            double maxDispl_, double maxTimeStep_,
                            NeighbourDetector<Agent> & neighbourDetector_,
                            std::ostream & logStream_=std::cerr)
        : unaryForces(unaryForces_), binaryForces(binaryForces_), bondForces(bondForces_),
          L(L_), maxDispl(maxDispl_), maxTimeStep(maxTimeStep_),
          neighbourDetector(neighbourDetector_),
          logStream(logStream_),
          totalIter(0), totalTime(0.0) { }

    /// Integrate the overdamped equation of motion for the system of agents in agentContainer.
    /// The stopping criterion is supplied by the IntegratorController object.
    virtual void integrate(std::vector<Agent*> & agents, IntegratorController & controller);

    int getTotalIter() const { return totalIter; }
    double getTotalTime() const { return totalTime; }

private:
    /// Advance solution by a time step no larger than dt;
    /// the time step is adapted to be consistent with maxDispl and maxTimeStep constraints.
    void advanceInitialStep(std::vector<Agent*> & agents, double * dt, double * displMax, util::SimpleStatistics & residualStatistics);
    void advanceSecondStep(std::vector<Agent*> & agents, double * dt, double * displMax, util::SimpleStatistics & residualStatistics);
    void advanceThirdStep(std::vector<Agent*> & agents, double * dt, double * displMax, util::SimpleStatistics & residualStatistics);
    void advanceFinalStep(std::vector<Agent*> & agents, double * dt, double * displMax, util::SimpleStatistics & residualStatistics);

    /// force and statistics calculation without updating positions
    void calculateForcesAndStress(std::vector<Agent*> & agents);
    void clearForces(std::vector<Agent*> & agents);

    void calculatePairwiseForce (std::vector<Agent*> & agents, size_t iAgent, Agent* agent2, std::vector<point_t> & force);
    void log(int totalIter, double totalTime, int iter, double time, double dt, util::SimpleStatistics const& residualStatistics);

private:
    UnaryForce <Agent> & unaryForces;
    BinaryForce <Agent> & binaryForces;
    BondForce <Agent> & bondForces;
    double L; //characteristic length
    double maxDispl, maxTimeStep;
    NeighbourDetector<Agent> & neighbourDetector;
    std::ostream & logStream;
    int totalIter;      // total number of iterations (all calls, needed for logging)
    double totalTime;   // total integration time (all calls, needed for logging)
    std::vector<point_t> force1, force2, force3, force4;
};

} // end namespace absmc

#endif
