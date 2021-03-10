#ifndef ABSMC_RUNGE_KUTTA_INTEGRATOR_HH
#define ABSMC_RUNGE_KUTTA_INTEGRATOR_HH

#include <algorithm>
#include <functional>
#include <limits>
#include <math.h>

#include "core/adaptiveRKIntegrator.h"
#include <omp.h>
#include <iostream>


namespace absmc {

    /// Requires force values in REAL UNITS (N)
    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::integrate(std::vector<Agent*> & agents, IntegratorController & controller)
    {
        double time = 0.0;
        int iter = 0;

        // make sure neighbourhood records are up to date
        neighbourDetector.updateSlowTimeScale(agents);
        // keep track of maximum displacement travelled by any agent at this step
        double displMax = 0.0;

        // statistics on force residuals (for convergence monitoring)
        util::SimpleStatistics residualStatistics;

        // initial guess for time step
        double dt = maxTimeStep;
        const size_t nAgents = agents.size();

        force1.resize(nAgents);
        force2.resize(nAgents);
        force3.resize(nAgents);
        force4.resize(nAgents);
        clearForces(agents);
        while (!controller.stop(iter, time, dt, residualStatistics)) {
            residualStatistics.reset();
            advanceInitialStep (agents, &dt, &displMax, residualStatistics);
            advanceSecondStep (agents, &dt, &displMax, residualStatistics);
            advanceThirdStep (agents, &dt, &displMax, residualStatistics);
            advanceFinalStep (agents, &dt, &displMax, residualStatistics);

            neighbourDetector.updateFastTimeScale(agents, displMax);
            time += dt;
            totalTime += dt;
            iter++;
            totalIter++;
            log(totalIter, totalTime, iter, time, dt, residualStatistics);
            //std::cout << "Step done, max force is " << residualStatistics.getMax() << std::endl;
        }
        calculateForcesAndStress(agents);
        //std::cout << "Integration converged" << std::endl;
    }

    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::log(int totalIter, double totalTime,
                                                                        int iter, double time,
                                                                        double dt,
                                                                        util::SimpleStatistics const& residualStatistics)
    {
        #ifdef ABSMC_INTEGRATOR_LOGGING
        logStream << totalIter << "  " << totalTime << " " << dt << " ";
        if (iter>0) {
            logStream << residualStatistics.getMin() << " "
                      << residualStatistics.getMax() << " "
                      << residualStatistics.getAvg() << " "
                      << residualStatistics.getNumValues() << " "
                      << "0 0 0 0 0 "
                      << residualStatistics.getL2Norm();
        }
        logStream   << std::endl;
        #endif
        //std::cout << "Total iterations " << totalIter << ", total time " << totalTime << ", dt " << dt << ", Residual max " << residualStatistics.getMax() << std::endl;
        //logger::info("Total iterations %d, total time %f, dt %f, Residual max %f", totalIter, totalTime, dt, residualStatistics.getMax());
    }

    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::advanceInitialStep(std::vector<Agent*> & agents,
                                                                                         double * dt, double * displMax,
                                                                                         util::SimpleStatistics & residualStatistics)
    {
        // maximum force exerted on any agent at this time step
        double forceMaxSqr = 0.0;

        const size_t nAgents = agents.size();
        point_t prevForce;
        double forceNormSqr;

        //calculate pairwise forces
        #pragma omp parallel for reduction(max : forceMaxSqr) private(prevForce, forceNormSqr) schedule(dynamic, 100)
        for (size_t iAgent=0; iAgent < nAgents; iAgent++) {
            Agent* agent = agents[iAgent];
            // skip obstacles
            if (!agent->isMobile) continue;
            const point_t prevForce = force1[iAgent];

            force1[iAgent] = unaryForces.calculateForce(*agent);

            if (agent->isAffected) {
                const size_t nNb = agent->neighbours.size();
                for (size_t iNb=0; iNb<nNb; iNb++) {
                    force1[iAgent] += binaryForces.calculateForce(*agent, *(agent->neighbours[iNb]) ) ;
                }

                const size_t nBonds = agent->bonds.size();
                for (size_t iNb=0; iNb<nBonds; iNb++) {
                    force1[iAgent] += bondForces.calculateForce(*agent, *(agent->bonds[iNb]) ) ;
                }
            }

            // keep track of maximum of forces exerted on agents
            forceNormSqr = normSqr<Agent::nDim>(force1[iAgent]);
            if (forceNormSqr>forceMaxSqr) {
                forceMaxSqr = forceNormSqr;
            }

            // compute force residuals
            for (size_t iDim=0; iDim<Agent::nDim; iDim++) {
                const double resF = fabs(force1[iAgent][iDim]-prevForce[iDim]);
                residualStatistics.takeValue(resF);
            }
        }

        const double forceMax = sqrt(forceMaxSqr);

        // adapt time step such that displacement is limited by maxDispl
        double newDt = *dt;
        if (!std::isnan(forceMax) ) {
            newDt = maxDispl / forceMax;
        }
        // make sure the new time step is no larger than the proposed one
        newDt = std::min(newDt, maxTimeStep);

        (*dt) = newDt;

        // maximum displacement travelled by any agent in this time step (real dimensions)
        (*displMax) = newDt * forceMax;

        point_t displacement;
        // update positions (real dimensions)
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
            // obstacle agents are not moved since their mobility equals zero
            for (size_t iDim = 0; iDim < Agent::nDim; iDim++) {
                agents[iAgent]->pos[iDim] +=
                    (agents[iAgent]->mobility[iDim]) * newDt * force1[iAgent][iDim] * 0.5;
            }
        }
        //std::cout << forceMax << std::endl;

    }



    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::advanceSecondStep(std::vector<Agent*> & agents,
                                                                                     double * dt, double * displMax,
                                                                                     util::SimpleStatistics & residualStatistics)
    {
        const size_t nAgents = agents.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (size_t iAgent=0; iAgent < nAgents; iAgent++) {
            Agent* agent = agents[iAgent];

            // skip obstacles
            if (!agent->isMobile) continue;

            force2[iAgent] = unaryForces.calculateForce(*agent);
            if (!agent->isAffected) continue;

            const size_t nNb = agent->neighbours.size();
            for (size_t iNb=0; iNb<nNb; iNb++) {
                force2[iAgent] += binaryForces.calculateForce(*agent, *(agent->neighbours[iNb]) ) ;
            }

            const size_t nBonds = agent->bonds.size();
            for (size_t iNb=0; iNb<nBonds; iNb++) {
                force2[iAgent] += bondForces.calculateForce(*agent, *(agent->bonds[iNb]) ) ;
            }
        }

        double newDt = *dt;
        // update positions (real dimensions)
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
            // obstacle agents are not moved since their mobility equals zero
            for (size_t iDim = 0; iDim < Agent::nDim; iDim++) {
                agents[iAgent]->pos[iDim] +=
                    (agents[iAgent]->mobility[iDim]) * newDt * (force2[iAgent][iDim] * 0.5 - force1[iAgent][iDim] * 0.5);
            }
        }
    }



    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::advanceThirdStep(std::vector<Agent*> & agents,
                                                                    double * dt, double * displMax,
                                                                    util::SimpleStatistics & residualStatistics)
    {
        const size_t nAgents = agents.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (size_t iAgent=0; iAgent < nAgents; iAgent++) {

            Agent* agent = agents[iAgent];

            // skip obstacles
            if (!agent->isMobile) continue;

            force3[iAgent] = unaryForces.calculateForce(*agent);
            if (!agent->isAffected) continue;

            const size_t nNb = agent->neighbours.size();
            for (size_t iNb=0; iNb<nNb; iNb++) {
                force3[iAgent] += binaryForces.calculateForce(*agent, *(agent->neighbours[iNb]) ) ;
            }

            const size_t nBonds = agent->bonds.size();
            for (size_t iNb=0; iNb<nBonds; iNb++) {
                force3[iAgent] += bondForces.calculateForce(*agent, *(agent->bonds[iNb]) ) ;
            }
        }

        double newDt = *dt;
        // update positions (real dimensions)
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
            // obstacle agents are not moved since their mobility equals zero
            for (size_t iDim = 0; iDim < Agent::nDim; iDim++) {
                agents[iAgent]->pos[iDim] +=
                    (agents[iAgent]->mobility[iDim]) * newDt * (force3[iAgent][iDim] * 1.0 - force2[iAgent][iDim] * 0.5);
            }
        }
    }


    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::advanceFinalStep(std::vector<Agent*> & agents,
                                                                                   double * dt, double * displMax,
                                                                                   util::SimpleStatistics & residualStatistics)
    {
        // maximum force exerted on any agent at this time step

        const size_t nAgents = agents.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {

            Agent* agent = agents[iAgent];

            // skip obstacles
            if (!agent->isMobile) continue;

            force4[iAgent] = unaryForces.calculateForce(*agent);
            if (!agent->isAffected) continue;

            const size_t nNb = agent->neighbours.size();
            for (size_t iNb=0; iNb<nNb; iNb++) {
                force4[iAgent] += binaryForces.calculateForce(*agent, *(agent->neighbours[iNb]) ) ;
            }

            const size_t nBonds = agent->bonds.size();
            for (size_t iNb=0; iNb<nBonds; iNb++) {
                force4[iAgent] += bondForces.calculateForce(*agent, *(agent->bonds[iNb]) ) ;
            }
        }
        double newDt = *dt;

        for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
            for (size_t iDim = 0; iDim < Agent::nDim; iDim++) {
                agents[iAgent]->pos[iDim] -=
                    (agents[iAgent]->mobility[iDim]) * newDt * (force3[iAgent][iDim] * 1.0);
                agents[iAgent]->pos[iDim] +=
                    (agents[iAgent]->mobility[iDim]) * newDt *
                    (force1[iAgent][iDim]+2.0*force2[iAgent][iDim]+2.0*force3[iAgent][iDim]+force4[iAgent][iDim])/6.0;
            }
        }
    }

    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::calculateForcesAndStress(std::vector<Agent*> & agents) {
        const size_t nAgents = agents.size();
        EuclideanMetrics<Agent::nDim, double> metrics;

        #pragma omp parallel for schedule(dynamic, 100)
        for (size_t iAgent=0; iAgent < nAgents; iAgent++) {

            agents[iAgent]->resetForce();
            agents[iAgent]->resetStressMatrix();

            Agent* agent = agents[iAgent];

            // skip obstacles
            if (!agent->isMobile) continue;

            unaryForces.accumulate(*agent);
            if (!agent->isAffected) continue;

            const size_t nNb = agent->neighbours.size();
            for (size_t iNb=0; iNb<nNb; iNb++) {
                point_t forceVec = binaryForces.calculateForce(*agent, *(agent->neighbours[iNb]) ) ;

                point_t const& pos1 = agent->pos;
                point_t const& pos2 = (agent->neighbours[iNb])->pos;

                // distance (real dimensions)
                const double distance = metrics.dist(pos1, pos2);
                const double R1 = agent->r;
                const double R2 = (agent->neighbours[iNb])->r;
                const double force = binaryForces.calculateForce(R1, R2, distance);


                const double xproject = (pos2[0]-pos1[0]) / distance;
                const double yproject = (pos2[1]-pos1[1]) / distance;
                const double zproject = (pos2[2]-pos1[2]) / distance;

                const double l1 = R1 + 0.5*(distance - R1 - R2);
                const double f12radeq = l1*force; //virial

                agent->addStress(point_t(xproject * f12radeq,   //not needed on intermediate steps
                                         yproject * f12radeq,       //used in the necrosis rule
                                         zproject * f12radeq) , R1*R1*R1);
                const double f_times_l_over_v = 0.75 * f12radeq/(3.1415926 * R1*R1*R1); // minus to make stress matrix components (and pressure) positive when forces are repulsive
                agent->addStressMatrix(f_times_l_over_v * xproject * xproject, // virial stress matrix, eq. to macroscopic Cauchy stress
                                        f_times_l_over_v * yproject * xproject,                    // elements are f_x * delta-r_x etc.
                                        f_times_l_over_v * zproject * xproject,
                                        f_times_l_over_v * yproject * yproject,
                                        f_times_l_over_v * zproject * zproject,
                                        f_times_l_over_v * zproject * yproject); // 00, 10, 20, 11, 22, 21

                agent->addForce(forceVec);
            }


            const size_t nBonds = agent->bonds.size();
            for (size_t iNb=0; iNb<nBonds; iNb++) {
                point_t forceVec = bondForces.calculateForce(*agent, *(agent->bonds[iNb]) ) ;

                point_t const& pos1 = agent->pos;
                point_t const& pos2 = (agent->bonds[iNb])->pos;

                // distance (real dimensions)
                const double distance = metrics.dist(pos1, pos2);
                const double R1 = agent->r;
                const double R2 = (agent->bonds[iNb])->r;
                const double force = bondForces.calculateForce(R1, R2, distance);


                const double xproject = (pos2[0]-pos1[0]) / distance;
                const double yproject = (pos2[1]-pos1[1]) / distance;
                const double zproject = (pos2[2]-pos1[2]) / distance;

                const double l1 = R1 + 0.5*(distance - R1 - R2);
                const double f12radeq = l1*force; //virial

                agent->addStress(point_t(xproject * f12radeq,   //not needed on intermediate steps
                                         yproject * f12radeq,       //used in the necrosis rule
                                         zproject * f12radeq) , R1*R1*R1);
                const double f_times_l_over_v = 0.75 * f12radeq/(3.1415926 * R1*R1*R1); // minus to make stress matrix components (and pressure) positive when forces are repulsive
                agent->addStressMatrix(f_times_l_over_v * xproject * xproject, // virial stress matrix, eq. to macroscopic Cauchy stress
                                        f_times_l_over_v * yproject * xproject,                    // elements are f_x * delta-r_x etc.
                                        f_times_l_over_v * zproject * xproject,
                                        f_times_l_over_v * yproject * yproject,
                                        f_times_l_over_v * zproject * zproject,
                                        f_times_l_over_v * zproject * yproject); // 00, 10, 20, 11, 22, 21

                agent->addForce(forceVec);
            }
        }
    }

    template<class Agent>
    void RungeKuttaIntegratorMPI<Agent>::clearForces(std::vector<Agent*> & agents) {
        const size_t nAgents = agents.size();
        for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
            agents[iAgent]->resetForce();
            agents[iAgent]->resetStressMatrix();
        }
    }

} // end namespace absmc

#endif
