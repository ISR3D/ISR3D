#ifndef ABSMC_AGENT_BASE_H
#define ABSMC_AGENT_BASE_H

#include <cstddef>
#include <iostream>
#include <cassert>
#include <vector>
#include <set>
#include <limits>

#include <util/agentTypeId.h>
#include "core/geometry.h"

namespace absmc {

template<size_t nDim_>
class AgentBase {
    template<size_t nDim>
    friend class AgentContainer;

    template<size_t nDim>
    friend class NBReader;

    template<typename AgentType>
    friend class Integrator;

    template<typename AgentType>
    friend class RungeKuttaIntegratorMPI;

    friend class BinaryForce3D;
    friend class BondForce3D;

    friend class Bilinear3D;
    friend class Wmlc3D;

    friend class WssGradientForce3D;
    friend class NeohookeanLJ3D;
    friend class RandomUnaryForce3D;

    friend class GrowthFactorForce3D;
    friend class UnaryStentDisplacementForce;

    class Bond {
    public:
        AgentBase* other;  //pointer to the other agent
        double length0;    //strain-free bond length

        Bond (AgentBase* other_, double length0_)
            : other(other_), length0(length0_) {}

        // get the bond's length
        // awkward that you have to pass this agent's pos here;
        // otherwise a pointer to this agent is needed for each bond,
        // and that bloats the memory
        double getDistance(Point<nDim_, double> pos) {
            return norm<nDim>(other->pos - pos);
        }

    };


public:
    typedef Point<nDim_, double> point_t;
    typedef std::vector<AgentBase*> nb_vec_t;
    typedef std::vector<Bond> bond_vec_t;

    static const size_t nDim = nDim_;
    static inline size_t _id = 0;

    AgentBase(AgentTypeId typeId_, point_t pos_, point_t  pos0_, double r_, point_t mobility_)
        : typeId(typeId_), pos(pos_), pos0(pos0_), r(r_), strain(0.0), wssOsi(0.0), wssMax(0.0) {
        setMobility(mobility_); //not a pure set
        id = _id++;
    }

    static size_t getMaxId() {
        return _id;
    }

    virtual AgentBase* clone() const = 0;

    //neighbour cleanup on deletion taken care of in agentContainer
    virtual ~AgentBase() {
        unsetSurfaceNormal();
        //clearNeighbours();
    }

    AgentTypeId getTypeId() const { return typeId; }

    //bool getGhostCondition() { return ghost; }
    //void setGhostCondition(bool ghost_) { ghost = ghost_; }

    point_t const& getPos() const { return pos; }
    void setPos(point_t const& newPos) { pos = newPos; }
    void setCoord(double const& newCoord, size_t const axis) { pos[axis] = newCoord; }

    point_t const& getPos0() const { return pos0; }
    void setPos0(point_t const& newPos0) { pos0 = newPos0; }

    double getR() const { return r; }
    void setR(double newR) { r=newR; }
    void setStrain(double newStrain)  {
        strain = newStrain;
    }

    double getStrain() const {
        return strain;
    }

    point_t const& getMobility() const { return mobility; }

    void setMobility(point_t const& newMobility) {
        mobility = newMobility;
        isMobile = false;
        for (size_t iDim=0; iDim<nDim; iDim++) {
            if(mobility[iDim] != 0.0)
                isMobile = true;
        }
    }

    inline bool getMobile() {
        return isMobile;
    }

    // whether the agent is affected by other agents in force calculation
    inline void setAffected(bool state) {
        isAffected = state;
    }

    inline bool const& getAffected() const {
        return isAffected;
    }

    /// write surface normal for agents at the surface
    inline void setSurfaceNormal(point_t const& normal) {
        surfaceNormal = new point_t(normal);
    }

    inline point_t const* getSurfaceNormal() const {
        return surfaceNormal;
    }

    /// remove the normal record
    inline void unsetSurfaceNormal() {
        delete surfaceNormal;
        surfaceNormal = nullptr;
    }

    inline bool isSurface() const {
        return (surfaceNormal != nullptr);
    }

    point_t const& getForce() const { return force; }
//    point_t const& getStressVec() const { return stressvec; }
    void resetForce() {
        force = point_t();
        stress = point_t();
        //strain = 0.0;
    }
    void resetStressMatrix() {
        stressmatrix[0] = 0.0;
        stressmatrix[1] = 0.0;
        stressmatrix[2] = 0.0;
        stressmatrix[3] = 0.0;
        stressmatrix[4] = 0.0;
        stressmatrix[5] = 0.0;
    }
    void resetpressure() {
        pressure=0.0;
    }
    inline void addForce(point_t const& newForce) {
            force += newForce;
    }

    //3D version of stress; Rcube is used to calculate dV inside the function
    void addStress(point_t const& newVirial, double Rcube) {
        double check = 0.0;
        for (size_t iDim=0; iDim<nDim; iDim++) {
            check = fabs(0.75*newVirial[iDim]/(3.1415926 * Rcube));
            stress[iDim] += check;
        }
    }

    void addStressMatrix(double s00, double s10, double s20, double s11, double s22, double s21) { // 00, 10, 20, 11, 22, 21

        stressmatrix[0] += s00;
        stressmatrix[1] += s10;
        stressmatrix[2] += s20;
        stressmatrix[3] += s11;
        stressmatrix[4] += s22;
        stressmatrix[5] += s21;

    }

    void addPressure(double pres) {
        pressure += pres;
    }
    double getPressure(){ return pressure;}

    double const& getStressMatrixComp(int comp) const { return stressmatrix[comp]; }

    point_t  const& getStress() const { return stress; }
    double getStressComp(int dim) const { return stressmatrix[dim]; }

    nb_vec_t & getNeighbours() { return neighbours; }
    nb_vec_t const& getNeighbours() const { return neighbours; }

    bond_vec_t & getBonds() { return bonds; }
    bond_vec_t const& getBonds() const { return bonds; }

    double getWssOsi() const { return wssOsi; }
    void setWssOsi(double wssOsi_) { wssOsi = wssOsi_; }

    double getWssMax() const { return wssMax; }
    void setWssMax(double wssMax_) { wssMax = wssMax_; }

    virtual void execAgentRule(size_t iAgent, std::vector<AgentBase*> & agentsCreated, std::vector<AgentBase*> & agentsDeleted) = 0;

    virtual void print(std::ostream & ostr) {
        ostr << typeId << " \t";
        for (size_t iDim=0; iDim<nDim; iDim++) ostr << pos[iDim] << " \t";
        for (size_t iDim=0; iDim<nDim; iDim++) ostr << pos0[iDim] << " \t";
        ostr << r << " \t";
    }

    virtual double getScalarByName(std::string const& quantity) const {
        if (quantity=="typeId") return typeId;
        else if (quantity=="posX") return pos[0];
        else if (quantity=="posY") { assert(nDim>=2); return pos[1]; }
        else if (quantity=="posZ") { assert(nDim>=3); return pos[2]; }
        else if (quantity=="pos0X") return pos0[0];
        else if (quantity=="pos0Y") { assert(nDim>=2); return pos0[1]; }
        else if (quantity=="pos0Z") { assert(nDim>=3); return pos0[2]; }
        else if (quantity=="radius") return r;
        else if (quantity=="mobilityX") return mobility[0];
        else if (quantity=="mobilityY") { assert(nDim>=2); return mobility[1]; }
        else if (quantity=="mobilityZ") { assert(nDim>=3); return mobility[2]; }
        else if (quantity=="forceNorm") { return norm<nDim>(force); }
        else if (quantity=="forceX") return force[0];
        else if (quantity=="forceY") { assert(nDim>=2); return force[1]; }
        else if (quantity=="forceZ") { assert(nDim>=3); return force[2]; }
        else if (quantity=="stressNorm") { return norm<nDim>(stress); }
        else if ((quantity=="stressX") || (quantity=="stressXX")) return stressmatrix[0];
        else if ((quantity=="stressY") || (quantity=="stressYY")) return stressmatrix[3];
        else if ((quantity=="stressZ") || (quantity=="stressZZ")) { assert(nDim>=3); return stressmatrix[4]; }
        else if (quantity=="stressXY") return stressmatrix[1];
        else if (quantity=="stressXZ") { assert(nDim>=3); return stressmatrix[2]; }
        else if (quantity=="stressYZ") { assert(nDim>=3); return stressmatrix[5]; }
        else if (quantity=="pressure") { assert(nDim>=3); return (- (stressmatrix[0] + stressmatrix[3] + stressmatrix[4]) / 3); }
        else if (quantity=="numNeighbours") return neighbours.size();
        else if (quantity=="numBonds") return bonds.size();
        else if (quantity=="strain") return strain;
        else if (quantity=="id") return id;
        else if (quantity=="isSurface") return isSurface();
        else if (quantity=="wssMax") return wssMax;
        else {
            // std::cout << "Warning: unknown scalar requested: " << quantity << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    virtual point_t getVectorByName(std::string const& quantity) const {
        if (quantity=="pos") return pos;
        else if (quantity=="pos0") return pos0;
        else if (quantity=="mobility") return mobility;
        else if (quantity=="force") return force;
        else if (quantity=="stress") return stress;
        else return point_t();
    }

    /// This method removes the agent from all neighbour lists and clears its own list
    void clearNeighbours() {
        for (typename nb_vec_t::iterator it = neighbours.begin() ; it != neighbours.end(); ++it) {
            (*it)->removeNeighbour(this);
        }
        neighbours.clear();
    }

    /// This method removes the agent from all bond lists and clears its own list
    void clearBonds() {
        for (typename bond_vec_t::iterator it = bonds.begin() ; it != bonds.end(); ++it) {
            (it->other)->removeBond(this);
        }
        bonds.clear();
    }

    inline void addBond(AgentBase* other) {
        addBond(other, true);
    }

    /// break the bond between this agent and another one
    inline void breakBond(AgentBase* other) {
        other->removeBond(this);
        removeBond(other);
    }

    /// This methods adds the cells from the neighbour list to bonds, up to maxBonds total bonds, located no further than maxDistance
    void addNeighboursToBonds(bool check = true, size_t maxBonds = std::numeric_limits<size_t>::max(), double maxDistance = std::numeric_limits<double>::max()) {
        size_t numBonds = bonds.size();
        if (numBonds >= maxBonds) return;

        //put all neighbours in a set sorted by distance
        auto cmp = [this](AgentBase* a, AgentBase* b) {
            return norm<nDim>(a->pos - this->pos) < norm<nDim>(b->pos - this->pos);
        };
        std::set< AgentBase* , decltype(cmp)> toBeAdded(cmp);
        for (auto& neighPtr : neighbours) {
            bool notAdded = true;
            if(check) {
                for (auto& bond : bonds) {
                    if (neighPtr == bond.other) {
                        notAdded = false;
                        break;
                    }
                }
            }
            if(notAdded) {
                toBeAdded.insert(neighPtr);
            }
        }
        for (auto& candidate : toBeAdded) {
            if ((numBonds >= maxBonds) || (norm<nDim>(candidate->pos - this->pos) > maxDistance)) break;
            if (candidate->bonds.size() < maxBonds ) {
                addBond(candidate);
                numBonds++;
            }
        }
    }

    // copy neighbours while keeping neighbourhoods symmetric
    void copyNeighbours(AgentBase const& other) {
        neighbours = other.neighbours;

        for (auto& nbr : neighbours) {
            auto& backLinks = nbr->getNeighbours();
            backLinks.push_back(this);
        }
    }

    size_t getId() const { return id; }
    void setId(size_t i) { id = i; _id = std::max(_id, i+1);} //usually a Very Bad Idea to call this directly, since it can lead to duplicate IDs; used for reading agents from file
    //int getGhostIndex() { return ghostindex; }
    //void setGhostIndex(int i) { ghostindex = i; }
private:

    //remove neighbour by address. Note that calling this by itself breaks the symmetry of neighbour lists
    void removeNeighbour(AgentBase* index) {
        for (typename nb_vec_t::iterator it = neighbours.begin() ; it != neighbours.end(); ++it) {
            if (*it == index) {
                neighbours.erase(it);
                return;
            }
        }
    }

    //remove bond by address. Note that calling this by itself breaks the symmetry of bonds
    void removeBond(AgentBase* index) {
        for (typename bond_vec_t::iterator it = bonds.begin() ; it != bonds.end(); ++it) {
            if (it->other == index) {
                bonds.erase(it);
                return;
            }
        }
    }

    inline void addBond(AgentBase* other, bool symmetric) {
        if (symmetric) {
            other->addBond(this, false);
        }
        bonds.push_back(Bond(other, norm<nDim>(other->pos - this->pos)));
    }

    inline void addBond(AgentBase* other, double length0, bool symmetric) {
        if (symmetric) {
            other->addBond(this, length0, false);
        }
        bonds.push_back(Bond(other, length0));
    }


    const AgentTypeId   typeId;
    //bool              ghost;
    point_t             pos, pos0;
    double              r, strain, pressure;
    point_t             mobility;
    bool                isMobile, isAffected = true;
    point_t*            surfaceNormal = nullptr;
    point_t             force;
//    point_t            stressvec;
    point_t             stress;
    nb_vec_t            neighbours;
    bond_vec_t          bonds;
    double              stressmatrix [6]; // 00, 10, 20, 11, 22, 21
    size_t              id;
    //int ghostindex;
protected:
    double wssOsi;
    double wssMax;

};

} // end namespace absmc

#endif
