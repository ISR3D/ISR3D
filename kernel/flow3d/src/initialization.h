#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include <iostream>
#include <cmath>
#include "flowFileInput.h"

/// A functional, used to instantiate boundary nodes at vessel openings
/*template<typename T>
class OpeningShapeDomain : public plb::DomainFunctional3D {
public:
    OpeningShapeDomain(plint xCoord, plb::plint cy_, plb::plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY) const {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
    }
    virtual OpeningShapeDomain<T>* clone() const {
        return new OpeningShapeDomain<T>(*this);
    }
private:
    plb::plint cx;
    plb::plint cy;
    plb::plint radiusSqr;
};
*/

/// A functional, used to initialize a pressure boundary to constant density
template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY, plint iZ) const {
        return density;
    }
private:
    T density;
};


// Look for the min and max of the slice defined by centerX_
template<typename T>
void findCenter(const MultiScalarField3D<int>& charMask, T& centerX_, T& centerY_, T& centerZ_, plint& curMinY, plint& curMinZ, plint& curMaxY, plint& curMaxZ) {
    curMinY = charMask.getNy() - 1;
    curMinZ = charMask.getNz() - 1;
    curMaxY = 0;
    curMaxZ = 0;

    /// this is very inefficient, since we are asking for a single byte from another node a lot of times
    /// should only be used for initial setup
    size_t b0counter = 0;
    size_t b1counter = 0;
    size_t b2counter = 0;
    for (plint i = 0; i < charMask.getNy(); i++) {
        for (plint j = 0; j < charMask.getNz(); j++) {
            int b = charMask.get(centerX_, i, j);
            //if the node is fluid
            if (b == 0) {
                b0counter++;
                if (i <= curMinY) {
                    curMinY = i;
                }
                if (i >= curMaxY) {
                    curMaxY = i;
                }
                if (j <= curMinZ) {
                    curMinZ = j;
                }
                if (j >= curMaxZ) {
                    curMaxZ = j;
                }
            }
            if (b == 1) {
                b1counter++;
            }
            if (b == 2) {
                b2counter++;
            }
        }
    }
    centerY_ = (curMinY + curMaxY) / 2; //might be between nodes for even number of nodes
    centerZ_ = (curMinZ + curMaxZ) / 2; //if the vessel isn't at z=0, find rightmost and leftmost empty cells and average

    //logger::info("At lowest X, %d fluid cells, %d solid cells, %d surface cells, %d cells total.", b0counter, b1counter, b2counter, parameters.getNy() * parameters.getNz());
    //logger::info("Center at %f, %f.", centerY, centerZ);
}


/// Generates a Poiseulle flow profile for a tube in XY plane
/// (assuming its main axis has z = 0, and is in the middle of the computational domain)
/// alpha is the angle between the tube's axis and OX (in radians)
/// Positive alpha correspond to a tube oriented toward positive X and Y
/// alphas outside (-pi/2, pi/2) can cause strange results
/// Values close to +- pi/2 may and will result in ridiculously large sections.
/// The current version is used ONLY to set the inlet boundary velocity
template<typename T, template<typename U> class Descriptor>
class PoiseuilleVelocity3D{
public:
    PoiseuilleVelocity3D(IncomprFlowParam<T>& parameters_,
                        const MultiScalarField3D<int>& charMask_)
        : parameters(parameters_), charMask (charMask_)
    {
        /// If this procedure is used in the future to set flow in the whole vessel,
        /// the center point should be replaced with a center line.
        centerX = 0;
        centerXnb = 7; /// arbitrary step to get a better estimation of alpha
        findCenter<T>(charMask, centerX,centerY,centerZ, lowestXminY, lowestXminZ, lowestXmaxY, lowestXmaxZ);

        plint nbXminY, nbXminZ, nbXmaxY, nbXmaxZ;
        findCenter<T>(charMask, centerXnb,centerYnb,centerZnb,nbXminY, nbXminZ, nbXmaxY, nbXmaxZ);

        //angle between the inlet pipe and the OX axis (in radians)
        //positive value means the inlet is facing in positive Y direction
        phi = atan2(centerYnb - centerY, centerXnb - centerX);

        //angle for z direction
        theta = acos(0.) - atan2(centerZnb - centerZ, centerXnb - centerX);

        pcout << "Inlet angle: phi " << phi << " radians, theta " << theta << " radians" << std::endl;
    }

    // returns the poiseuille velocity; for creating boundary conditions, since the profile is only correct for X=0
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const {

    //calculate the velocity magnitude in (iX, iY, iZ)
        T uMagn;
        /// TODO: add theta to radius calculation
        const T r_sq = pow((iY - centerY)* cos(phi), 2) + pow((iZ - centerZ), 2);
        const T R_sq = (pow((lowestXmaxY - centerY)* cos(phi), 2) + pow((lowestXmaxZ - centerZ), 2)) / 2;
        const T ratio_sq = r_sq / R_sq;
        // may in principle set velocity for solid nodes
        if (ratio_sq < 1)// && (!b))
        {
            uMagn = parameters.getLatticeU() * (1 - ratio_sq);
            ////logger::info("Border velocity set to %g", uMagn);
            ////logger::info("Border velocity along X set to %g", uMagn * cos(alpha));
            ////logger::info("Border velocity along Y set to %g", uMagn * sin(alpha));
        } else
        {
            uMagn = 0;
            ////logger::info("Border velocity set to zero");
        }
        u[0] = uMagn * sin(theta) * cos(phi);
        u[1] = uMagn * sin(theta) * sin(phi);
        u[2] = uMagn * cos(theta);
    }

private:

    IncomprFlowParam<T> parameters;
    T phi, theta; // polar coords, ISO 80000-2:2019 naming
    T centerX, centerY, centerZ;
    T centerXnb, centerYnb, centerZnb; // center point of a neighbouring slice to determine the angle
    plint lowestXminY, lowestXminZ, lowestXmaxY, lowestXmaxZ; // we need channel size to compute flow profile

    const MultiScalarField3D<int>& charMask;
};


//returns Poiseulle velocity for a cylinder oriented along the X axis
//for a very rough initial approximation of velocity
//Also returns a linearly dropping density (same as pressure for incompressible LB) from inlet to outlet
template <typename T>
class CylindricalPoiseuilleDensityAndVelocity
{
public:

    CylindricalPoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_)
    : parameters(parameters_) { }

    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T, 3 > & u) const
    {
        u[0] = T();//poiseuilleVelocity(iY, iZ);
        u[1] = T();
        u[2] = T();

        const T Lx = parameters.getNx()-1;
        const T uMax = parameters.getLatticeU();
        const plint N = parameters.getResolution();
        //rho = (T) 1;
        //Linearly dropping rho for better initial numerical convergence
        const T poiseuillePressure = 8.* uMax*uMax / parameters.getRe()
                / N * (Lx/(T)2-(T)iX);
        /// Convert pressure to density according to ideal gas law
        rho = poiseuillePressure*DESCRIPTOR<T>::invCs2 + (T)1;

        //1.0 - DESCRIPTOR<T>::invCs2 * (dt*dt)/(dx*dx)* 48.0/parameters.getRe() * (iX-(nx-1)/2)/(nx-1);
        //1.0 + invCs2 * (1)/(v^2)* 48.0/parameters.getRe() * ((L)/2-iX)/(L);
    }
private:
    IncomprFlowParam<T> parameters;

    T poiseuilleVelocity(const plint iY, const plint iZ) const
    {

        const plint ny = parameters.getNy();
        const plint nz = parameters.getNz();

        //roughly equal to the vessel diameter
        const plint N = parameters.getResolution();

        const T r_sq = pow((iY - (T) (ny - 1) / 2), 2) + pow((iZ - (T) (nz - 1) / 2), 2);
        const T R_sq = (T)(N-2)*(T)(N-2)/4;//trying to adjust to real geometry taking into account # of fluid nodes
        const T y_sq = r_sq / R_sq;
        if (y_sq <= 1)
        {
            return parameters.getLatticeU() * (1 - y_sq);
        } else
        {
            return 0;
        }
    }

};

///This approximation is used as a boundary condition for rectangular channel
template <typename T>
class TwoPlatePoiseuilleVelocity {
public:
    TwoPlatePoiseuilleVelocity(IncomprFlowParam<T> parameters_,
                        const MultiScalarField3D<int>& charMask_)
        : parameters(parameters_), charMask (charMask_)
    {

        //find the center at lowest-X slice
        centerX = 0;
        centerZ = (parameters.getNz() - 1) / 2; //if the vessel isn't at z=0, find rightmost and leftmost empty cells and average
        top = parameters.getNy() - 1;
        bottom = 0;
        for (int i = 0; i <= parameters.getNy(); i++) {
            int b  = charMask.get(centerX, i, centerZ);
            ////logger::info("b = %d", b);
            if (b == 0) {
                //logger::info("Bottom border of inlet at node %d", i);
                bottom = i;
                break;
            }
        }
        //logger::info("%d nodes along axis Y", parameters.getNy());
        for (int i = parameters.getNy() - 2; i >= 0; i--) {
            int b  = charMask.get(centerX, i, centerZ);
            if (b == 0) {
                //logger::info("Top border of inlet at node %d", i);
                top = i;
                break;
            }
        }
        centerY = (top + bottom) / 2; //might be between nodes for even number of nodes
    }

    // returns the poiseuille velocity; for creating boundary conditions
    // rather slow, since it looks for the bounds each time
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const {

    //calculate the velocity magnitude in (iX, iY, iZ)
        T uMagn;
        //distance along Y
        const T r_sq = pow((iY - centerY), 2);
        const T R_sq = pow((top - centerY), 2);
        const T ratio_sq = r_sq / R_sq;
        if (ratio_sq < 1)
        {
            uMagn = parameters.getLatticeU() * (1 - ratio_sq);
            ////logger::info("Border velocity set to %g", uMagn);
            ////logger::info("Border velocity along X set to %g", uMagn * cos(alpha));
            ////logger::info("Border velocity along Y set to %g", uMagn * sin(alpha));
        } else
        {
            uMagn = 0;
            ////logger::info("Border velocity set to zero");
        }
        u[0] = uMagn;
        u[1] = T();
        u[2] = T();
    }

private:
    IncomprFlowParam<T> parameters;
    T centerX;
    T centerY;
    T centerZ;
    plint top, bottom;
    const MultiScalarField3D<int>& charMask;
};


#endif // INITIALISATION_H
