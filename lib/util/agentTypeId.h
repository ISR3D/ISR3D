#ifndef ABSMC_TYPE_ID_H
#define ABSMC_TYPE_ID_H

#include <iostream>

namespace absmc {

/// Make sure to keep the MAX element up to date when adding new types
enum AgentTypeId { tAny, tSMC2D, tIEL2D, tObstacle2D, tSMC3D, tIEL3D, tObstacle3D, tEEL3D, tEndo3D, tECM3D, AgentTypeId_MAX = tECM3D};

inline std::ostream & operator<<(std::ostream & ostr, const AgentTypeId typeId)
{
    return ostr << static_cast<int>(typeId);
}

inline std::istream & operator>>(std::istream & istr, AgentTypeId& typeId)
{
    int typeIdInt;
    istr >> typeIdInt;
    typeId = static_cast<AgentTypeId>(typeIdInt);
    return istr;
}


} // end namespace absmc

#endif
