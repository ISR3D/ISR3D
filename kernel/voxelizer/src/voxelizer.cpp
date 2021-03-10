#include <cstdlib>
#include <list>
#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>
#include "voxelizer.hpp"

using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

void voxelizer(int argc, char * argv[]) {

    Instance instance(argc, argv, {
            {Operator::F_INIT, {"cellPositionsIn"}},
            {Operator::O_F, {"domainOut1",
                             "domainOut2",
                             "cellsOut",
                             "surfacePositionOut",
                             "domainBoxOut"}}
    });

    VoxelizingMapper voxelizingMapper(instance);

    while (instance.reuse_instance()) {
        voxelizingMapper.execute();
    }
}


int main(int argc, char * argv[]) {
    voxelizer(argc, argv);
    return EXIT_SUCCESS;
}
