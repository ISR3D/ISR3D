#ifndef ENDO_RULES_3D_H
#define ENDO_RULES_3D_H

#include <vector>
#include <cassert>
#include <iostream>

#include "core/agentRule.h"
#include "model3D/endo3D.h"
#include "core/geometry.h"
#include <util/agentTypeId.h>

namespace absmc {


	/// Endo cell cycle rule. This is SMC growth rule with brakes (NO inhibition) removed
	class EndoCellCycleRule3D : public AgentRule<Endo3D> {
	public:
		typedef Endo3D agent_t;
		typedef AgentBase<3> agent_base_t;
		typedef Point<3, double> point_t;

		EndoCellCycleRule3D(double drugConcThreshold_, double dt_)
        : drugConcThreshold(drugConcThreshold_), dt(dt_) { }



		virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

			const int tCheckPoint = 0;  // G0 checkpoint (at beginning of cell cycle)
			const int tMitosis = agent.lengthG1 + Endo3D::lengthSG2M; //hours
			const double t = agent.getClock(); //hours
			const double growthFactor = pow(2.0, 1.0/(3.0*agent.lengthG1/dt));

			// at G0 check point or in G1
			if (t < agent.lengthG1) {
				// check for: biological activity, drug concentration, wss threshold, contact inhibition.
				// if any criterion for prohibition is met, G0.
				if (!agent.getBaFlag() || agent.getDrugConc() >= drugConcThreshold || agent.getCiFlag() || agent.getNOFlag())
				{
					//return to G0 and stop proliferation
					if (t > 0) {
						agent.setState(Endo3D::G0);
						agent.setClock(tCheckPoint);
						agent.setBaFlag(false);
					}
				}
				// otherwise enter G1 and increase clock
				else {
					agent.setState(Endo3D::G1);
					agent.setClock(t+dt);
					agent.grow(growthFactor);
				}
			}
			// S/G2/M
			else {
				if (t < tMitosis) { // cell in S/G2/M
					agent.setState(Endo3D::SG2M);
					agent.setClock(t+dt);
				}
				else { // cell at M: perform mitosis

					// random displacement from position of parent cell
					// TODO:stay on the surface while displacing. Needs a surf. normal for that?
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
					agent.setState(Endo3D::G0);
					agent.setClock(0);
					agent.setAge(agent.getAge()+1);
					agent.getNeighbours().clear();


					Endo3D* agentB = new Endo3D(posB, posB, newR);
					agentB->setAgentRule(agent.getAgentRule() );
					agentB->getNeighbours().clear();

					agentsCreated.push_back(agentB);

				}
			}
		}

	private:
		const double drugConcThreshold;
		const double dt;
	};


	class EndoNecrosisRule3D : public AgentRule<Endo3D> {
	public:
		typedef Endo3D agent_t;
		typedef AgentBase<3> agent_base_t;
		typedef Point<3, double> point_t;
		typedef std::vector<agent_base_t*> nb_vec_t;
		typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

		EndoNecrosisRule3D(double maxStress_, double maxR) : maxStress(maxStress_), maxRsqr(maxR * maxR) { }

		virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {

			//if the agent isn't attached to anything, delete it
			nb_vec_t nbs = agent.getNeighbours();
			const size_t nNb = nbs.size();
			size_t nClose = 0;
			for (size_t iNb=0; iNb<nNb; iNb++) {
				if(metrics.distSqr(agent.getPos(), nbs[iNb]->getPos()) <= maxRsqr) {
					nClose++;
				}
			}
			if (nClose == 0) {
				std::cout << "Endothelium necrosis rule, deleting agent: no neighbours" << std::endl;
				agentsDeleted.push_back(&agent);
				return;
			}
			// if stress is larger than threshold, delete agent
			if (norm<3>(agent.getStress() ) > maxStress) {
				std::cout << "Endothelium necrosis rule, deleting agent: stress=" << norm<3>(agent.getStress()) << " > maxStress" << std::endl;
				agentsDeleted.push_back(&agent);
				return;
			}
		}

	private:
		double maxStress;
		double maxRsqr;
		metrics_t metrics;
	};


} // end namespace absmc

#endif
