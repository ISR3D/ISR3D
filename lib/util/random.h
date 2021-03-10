#ifndef ABSMC_RANDOM_H
#define ABSMC_RANDOM_H

#include <cstdlib>
#include <cmath>

namespace util  {

class Random {
public:
    static void init(unsigned int seed)
    {
        srand(seed);
    }
    static void init()
    {
        srand(time(NULL));
    }
    /// return uniformly distributed floating-point pseudo-random number from range [0,max]
    static double getNext(double max) { return max * rand() / (double)(RAND_MAX); }

    /// return pseudo-Gaussian distributed pseudo-random number with mean x0 and width sigma
    // all values lie within 3 sigma from x0
    static double getNextGauss(double x0, double sigma) {
        // straightforward acceptance-rejectance method, very inefficient
        const double s3 = 3.0*sigma;
        const double c = -0.5/(sigma*sigma);
        double x,y;
        do {
            x = s3*(2*getNext(1.0)-1);
            y = getNext(1.0);
        } while (y > exp(c*x*x) );
        return x+x0;
    }
};

} // end namespace util

#endif
