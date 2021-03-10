#ifndef ABSMC_MATH_UTIL_H
#define ABSMC_MATH_UTIL_H

#include <cmath>

namespace absmc {

namespace math  {

// continental pi
const double pi = acos(-1.0);

inline double sgn(const double x)
{
    return (x >= 0) ? 1.0 : -1.0;
}

} // end namespace math

} // end namespace absmc

#endif
