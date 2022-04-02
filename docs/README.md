![Diagram](https://github.com/ISR3D/ISR3D/blob/master/docs/SchematicDiagram.jpg =100x20)

Contents
----------

* [About ISR3D & MUSCLE3](#About-ISR3D-&-MUSCLE3)
* [Model description](#Model-description)
* [Quick installation guide](#Quick-installation-guide)
* [Run the simulation](#Run-the-simulation)
* [Uncertainty quantification](#Uncertainty-quantification)
* [Visualization](#Visualization)
* [Uniaxial strain test](#Uniaxial-strain-test)


## About ISR3D & MUSCLE3
The «3D in-stent restenosis» (ISR3D) model is released as part of the InSilc project (https://insilc.eu). Currently two versions have been released. Each of them is associated to a corresponding publication. You can find the exact version of the ISR3D used in that pulication in our [release](https://github.com/ISR3D/ISR3D/releases).

ISR3D is a multiscale model simulating the post-stenting tissue growth. It consists of multiple submodels and the communication between the submodels is build by MUSCLE (Multiscale Coupling Library and Environment). We invite user to [MUSCLE3](https://muscle3.readthedocs.io/en/latest/) for more details. 

## Model description
The ISR3D application is intended for simulation of smooth muscle cells (SMC) proliferation and restenosis process as a complication
after coronary stenting procedure. The 3D in-stent restenosis model (ISR3D) is a fully coupled 3D multiscale
model, which includes several single-scale submodels as well as utility modules which facilitate communication
between the submodels. The submodel structure is similar to the one used by Caiazzo et al. in [1] and is described in detail in [2–3]. The submodels are described in this section.

**Mechanical model of SMCs.** The model code is located in ``kernel/absmc``. The agent-based model is used for modelling of mechanical response of vessel’s
walls. The cells of vessel tissue and stent struts are presented as an array of point-wise particles with no mass, and
the interactions between them are provided by potential forces of repulsion and adhesion. The effective radii of
particles represent the radii of corresponding cells and change during growth. Mathematically the problem is
formulated as a Cauchy problem for a system of ordinary differential equation. The coordinates of agents after each
step of biological model or stent deployment model are taken as initial conditions.

**Biological model of SMCs** is located in ``kernel/absmc`` as well. Cell cycle model is used for modelling the cell dynamics. Cell lifecycle is a sequence of
growth, replication and division of the cell; at the end of the lifecycle the cell divides into two daughter cells. The
processes that influence the cell lifecycle take place in the 30 μm neighbourhood around the cell; time scale of one
cycle is around 24-48 hours.

The growth of separate cells is modelled by a finite-state automaton. For each cell it can be in the state of growth
(G1), synthesis/repeat growth/mitosis (S/G2/M), or be idle (G0). Cells move from one state to the next, and stop or
die under the influence of external factors such as mechanical stresses (from mechanic model of SMC), the
concentration of nitric oxide (calculated from the shear stresses passed by the hydrodynamic model), and the
concentration of growth suppressing drugs (from diffusion model). The biological model provides new radii, states
and coordinates of the cells as its output. Growth of the neointima takes several dozens of cell cycles and stops
several weeks after the stenting procedure. Growth can progress up to the full vessel occlusion in some special
cases.

**Hydrodynamic model.**  The model code is located in ``kernel/flow3d``. Blood flow in the stented vessel is driven by a stationary Newton hydrodynamic model,
which provides relevant range of shear stresses on the vessel walls. The simulation of blood flow in the vessel is
provided by the lattice Boltzmann calculation method in 3D rectangular mesh (D3Q18). The field of flow velocities in each
lattice cell at the next timestep, and shear stress values in the cells near the wall of the vessel are the output of the
solver.

**Utility mapping modules.** To produce the input data for hydrodynamic solver the array of agents has to be
presented as a surface of the vessel wall. Each cell of hydrodynamic solver mesh where agent is present is marked as
solid, and then the obtained configuration is smoothed, which is done by the ``voxelizer`` module. The obtained voxel domain is then passed to the ``flow3d`` solver (and potentially other models acting on a 3D grid) by ``distributor``. The values of shear stresses are mapped back to the corresponding agent-cells by the ``collector``.

[1] Caiazzo, A., Evans, D., Falcone, J. L., Hegewald, J., Lorenz, E., Stahl, B., Wang, D., Bernsdorf, J., Chopard, B., Gunn, J. P., Hose, D. R., Krafczyk, M., Lawford, P. V., Smallwood, R., Walker, D., & Hoekstra, A. G. (2011). A Complex Automata approach for in-stent restenosis: Two-dimensional multiscale modelling and simulations. Journal of Computational Science, 2(1), 9–17. https://doi.org/10.1016/j.jocs.2010.09.002

[2] Zun, P. S., Anikina, T., Svitenkov, A., & Hoekstra, A. G. (2017). A Comparison of Fully-Coupled 3D In-Stent Restenosis Simulations to In-vivo Data. Frontiers in Physiology, 8(May), 284. https://doi.org/10.3389/fphys.2017.00284

[3] Zun, P. S., Narracott, A. J., Chiastra, C., Gunn, J., & Hoekstra, A. G. (2019). Location-Specific Comparison Between a 3D In-Stent Restenosis Model and Micro-CT and Histology Data from Porcine In Vivo Experiments. Cardiovascular Engineering and Technology, 10(4), 568–582. https://doi.org/10.1007/s13239-019-00431-4

For additional information please contact Pavel Zun, <pavel.zun (at) gmail.com>


## Quick Installation Guide

### MUSCLE3

The source code can be obtained from https://github.com/multiscale/muscle3/

At this point, it is recommended to use MUSCLE3 version 0.4.0 available at https://github.com/multiscale/muscle3/releases/tag/0.4.0.

The installation process is described in detail in the MUSCLE3 guide, located at https://muscle3.readthedocs.io/en/latest/installing.html

Both Python and C++ parts of MUSCLE3 are required. The C++ part should be built with MPI support to enable parallel flow computation.

NB: before running models with MUSCLE3 the library path has be loaded by running

``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<MUSCLE3 installation location>/lib``

(This can also be added to your ``~/.bashrc`` to run automatically on terminal startup)

If you have installed MUSCLE3 in a Python virtual environment, as the guide recommends, the environment also has to be activated by running 

``source <MUSCLE3 venv location>/bin/activate``

### ISR3D

Dependencies: in addition to Muscle, ``cmake`` ver. 3.6.3 or later. Additional requirements are a MPI-compatible C++14 compiler (e.g. GCC 6+ and OpenMPI).

The ISR3D model is capable of working with different flow solvers. This version is set up to use Palabos library, developed by UNIGE. By default, if Palabos is not detected when building the flow solver, it will be downloaded from the UNIGE repository using the script located at ``ISR3D/lib/install_palabos.sh``, which will install Palabos 2.2.0 to the ``/lib/`` folder.

To compile ISR3D and install it into the active directory, create and run ``./ISR3D/build.<machinename>.sh``
An example of this file for a desktop Debian-like machine is ``build.linux.sh``. Make sure to set ``MUSCLE3_HOME`` to the location of your MUSCLE3 installation.
For different machines, especially in HPC environments, it is necessary to prescribe the appropriate compiler name and the address of MUSCLE3 (and flags, if necessary). The folder also contains an example for the Dutch national supercomputer Cartesius, ``build.cartesius.sh`` and the Polish supercomputer Eagle, ``build.eagle.sh``.

After running the build script, you should see a CMake build output for each submodel, ending with ``BUILD SUCCESSFUL``, and the executables can be then be found in ``ISR3D/build``.
To clean the build folders and remove CMake cache files, run ``./build.<machinename>.sh clean``. Cleaning the cache is recommended whenever your library configuration changes.


## Run the simulation
To test MUSCLE3 (on its own), please follow the guide for running test examples located at https://muscle3.readthedocs.io/en/latest/cplusplus.html.

The simulation of the ISR3D can be divided into two parts. The first part is about initial stent deployment. The stent is expanded radially with a capsule-shaped balloon until it reaches a predefined deployment depth. The parameters needed for the process are all listed in a [configuration file](https://github.com/ISR3D/ISR3D/tree/master/config). To run this part, you can use following commands:
```
cd kernel/absmc
./prepareGeneratedStent.sh /directory-to-cfg-file/input_stage1.cfg
```

At the end of the process, a multiple stage3 files will be generated. There are files are the necessary inputs for the second part of the simulation. The second part is about post-stenting tissue growth. The parameters of this process are all listed in [ymmsl configuration file](https://github.com/ISR3D/ISR3D/tree/master/config). To run the simulation, open several separate terminals, or through &, or use a bash script, run:
```
muscle_manager /directory-to-stage4-ymmsl-file/input_stage4.ymmsl
/Installation-directory-of-ISR3D/build/smc --muscle-instance=smc
/Installation-directory-of-ISR3D/build/voxelizer --muscle-instance=voxelizer
/Installation-directory-of-ISR3D/build/distributor --muscle-instance=distributor
mpirun -np 2 /Installation-directory-of-ISR3D/build/FlowController3D --muscle-instance=flow
/Installation-directory-of-ISR3D/build/collector --muscle-instance=collector
```
Don't forget to load the MUSCLE3 library and enable venv if you use it. If you have done everything correctly, the models should start producing output as soon as all of them connect to the manager.

## Uncertainty quantification
To perform uncertainty quantification, a large number of evaluations are needed. ISR3D are fairly computational expensive. With a small vessel (2mm diameter), it takes around 500 to 600 core-hour for a single evaluation depending on the amount of growth. You would need some supercomputer resource to carry it out. We provide a [python script](https://github.com/ISR3D/ISR3D/blob/master/scripts/sample_generator.py) by which samples (Sobol sequence) of a UQ compaign can be generated and write them into corresponding configuration file (.cfg and .ymmsl). 

## Visualization
The .vtp output files can be visualized with Paraview. Please note that by default Paraview does not show point data, so after loading the file you have to use the "Select points through" option and select the data area to visualize the points. After that you can use "Extract selection" to further manipulate the selected points.

It is also possible to process the output result with VMTK.


## Uniaxial strain test
For uninaxial strain tests, first a Python script is used to generate the tissue sample, then it is converted to a .dat format, equilibrated, and finally the uniaxial strain test itself is performed. Naturally, the first three steps (up to equilibration) can be shared between multiple strain tests, as long as the forces' equilibrium distance is the same.

**Installation note:** since this test only involves the agent-based model, MUSCLE3 is not actually used. The necessary files can be built without relying on MUSCLE3. For this, in ``kernel/absmc/build.sh`` change the parameter to ``-DBUILD_MULTISCALE:BOOL=FALSE`` and run the ``./build.sh compile  && ./build.sh install`` directly. In this case, you still have to set up your compiler and the exports which are normally set in the machine-specific ``build.<machinename>.sh`` and the top-level ``ISR3D/build.sh``.

First, the sample has to be generated with a Python script. For this, Poisson disk sampling is used. The script for generation is located in ``scripts/`` and takes the output file address and the dimensions along X, Y, and Z axes (in mm) to generate a rectangular sample. For example:
```
cd scripts
python3 generate-cuboid-poisson-disc.py ../cxa/strain_test_sample.csv 2.0 0.4 1.2
```

The simulation parameters are denoted in a .cfg file. An example is located at ``cxa/stage1.strain_test.cfg``. The most notable parameters are the force coefficients ``strain_force_c1..c6``, the convergence level ``strain_convergence_level``, the maximum strain to be tested ``strain_max_value`` and the step between the strain values being tested ``strain_deform_step``.

Once the .csv file is ready, you have to convert it to .dat to use in the simulation. For this, 
```
cd ISR3D/kernel/absmc/build
./csvArteryLoader3D ../../../cxa/stage1.strain_test.cfg
./equilibration3D ../../../cxa/stage1.strain_test.cfg
```
After this, you can perform the uniaxial strain tests with 
```
./uniaxialStrainTest ../../../cxa/stage1.strain_test.cfg
```
The stress values for each strain value will be printed to the console.
