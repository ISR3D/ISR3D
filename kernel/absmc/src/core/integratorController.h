#ifndef ABSMC_INTEGRATOR_CONTROLLER_H
#define ABSMC_INTEGRATOR_CONTROLLER_H

#include <util/simpleStatistics.h>
#include <math.h>

namespace absmc {

/// Base class for integrator controllers.
/// Controls ode integrators by deciding when integration should be stopped.
class IntegratorController {
public:
    virtual ~IntegratorController() { }

    virtual bool stop(int iter, double time, double dt, util::SimpleStatistics const& residualStatistics) = 0;
};


/// Integration over an interval of fixed length.
class FixedIntervalController : public IntegratorController {
public:
    FixedIntervalController(double maxTime_) : maxTime(maxTime_) { }

    virtual bool stop(int iter, double time, double dt, util::SimpleStatistics const& residualStatistics) {
        if (time >= maxTime) {
//             std::cerr << time << "  " << iter << "  " << 1 << std::endl;
            return true;
        }
        return false;
    }
private:
    const double maxTime;
};

/// Integration over a fixed number of iterations.
class FixedIterationNumController : public IntegratorController {
public:
    FixedIterationNumController(double maxIter_) : maxIter(maxIter_) { }

    virtual bool stop(int iter, double time, double dt, util::SimpleStatistics const& residualStatistics) {
        if (iter >= maxIter) {
//             std::cerr << time << "  " << iter << "  " << 1 << std::endl;
            return true;
        }
        return false;
    }
private:
    const double maxIter;
};

/// Integration until maximum norm of the force residual vector is (instantaneously) below a defined threshold,
/// or time exceeds maximum integration time.
class ForceResidualMaxNormController : public IntegratorController {
public:
    ForceResidualMaxNormController(double maxTime_, double charForce_, double eps_)
        : maxTime(maxTime_), charForce(std::abs(charForce_)), eps(eps_) { }

    virtual bool stop(int iter, double time, double dt, util::SimpleStatistics const& residualStatistics) {
        if (time >= maxTime) {
//             std::cerr << time << "  " << iter << "  " << 1 << std::endl;
            return true;
        }
        if (iter>0 && residualStatistics.getMax() < eps*charForce) {
//             std::cerr << time << "  " << iter << "  " << 2 << std::endl;
            return true;
        }
        return false;
    };
private:
    const double maxTime;
    const double charForce; // system-characteristic force used to normalize residuals
    const double eps;       // desired convergence level
};


} // end namespace absmc

#endif
