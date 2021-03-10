#include <iostream>
#include <cmath>
#include <stdexcept>
#include "palabos3D.h"
#include "palabos3D.hh"
#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>
#include <util/logger.h>
#include "FlowSolver.h"

using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

int main(int argc, char * argv[]) {
    try
    {
        plbInit(&argc, &argv);
        bool plbInitialized = false;

        Instance instance(argc, argv, {
                {Operator::F_INIT, {
                    "geometryIn"}},    // grid of int
                {Operator::O_F, {
                    "shearStressOut"}} // grid of double
            },
            plb::global::mpi().getGlobalCommunicator(), plb::global::mpi().bossId()
        );

        FlowSolver flowSolver(instance);
        pcout << "Initializing flow solver..." << std::endl;

        cout << "Flow solver initialized at rank " << plb::global::mpi().getRank()  << std::endl;

        plb::global::mpi().barrier();

        pcout << "All flow solver processes initialized" << std::endl;


        while (instance.reuse_instance()) {
            pcout << "Receiving geometry" << std::endl;
            plb::global::mpi().barrier();
            flowSolver.receiveMask();
            plb::global::mpi().barrier();

            if (!plbInitialized) {
                /// initialize the flow model
                flowSolver.initialize();
                plbInitialized = true;
            }

            flowSolver.execute(argc, argv);
        }
    } catch (std::runtime_error& e)
    {
        logger() << "RUNTIME ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...)
    {
        logger() << "UNKNOWN ERROR" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

