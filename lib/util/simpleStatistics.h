#ifndef ABSMC_SIMPLE_STATISTICS_H
#define ABSMC_SIMPLE_STATISTICS_H

#include <math.h>
#include <limits>
#include <algorithm>

namespace util {

/// Compute statistics for a set of scalar values, updating statistics quantities instantaneously as values are added.
/// The number of values is not limited (except by the maximal value of int).
class SimpleStatistics {
public:
    SimpleStatistics() : totalCount(0) { reset(); }

    void reset() {
        minY =  std::numeric_limits<double>::max();
        maxY = -std::numeric_limits<double>::max(),
        sumY = 0.0;
        sumSqrY = 0.0;
//         sumAbsY = 0.0;
        numValues = 0;
    }

    void takeValue(double y) {
//        minY = std::min(minY, y);
        maxY = std::max(maxY, y);
//        sumY += y;
//        sumSqrY += y*y;
//         sumAbsY += fabs(y);
//        numValues++;
//        totalCount++;
    }

    double getMin() const { return minY; }
    double getMax() const { return maxY; }
    double getSum() const { return sumY; }
    double getL2Norm() const { return (sqrt(sumSqrY) / numValues); }
    double getAvg() const { return sumY / numValues; }
//     double getAbsSum() const { return sumAbsY; }
//     double getAbsAvg() const { return sumAbsY / numValues; }
    int getNumValues() const { return numValues; }
    int getTotalCount() const { return totalCount; }

private:
    double minY, maxY, sumY, sumSqrY;
//     double  sumAbsY;
    int numValues;
    int totalCount;
};

} // end namespace util

#endif
