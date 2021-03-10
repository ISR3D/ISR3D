#ifndef ABSMC_VTK_WRITER_H
#define ABSMC_VTK_WRITER_H

#include <string>
#include <vector>

#include "core/geometry.h"

namespace absmc {

namespace graphics {

template<class Agent>
class VtkWriter {
public:
    typedef Agent agent_t;
    typedef Point<Agent::nDim, double> point_t;

    static void writeVtkPolyData(std::vector<agent_t*> const& agents,
                                 std::string const& fileName,
                                 std::string const& scalarQuantityString="",
                                 std::string const& vectorQuantityString="");

private:
    static double convertNaNRep(double scalar, double NaNRep=0.0);
    static point_t convertNaNRep(point_t const& vector, double NaNRep=0.0);
};

} // end namespace graphics

} // end namespace absmc

#endif
