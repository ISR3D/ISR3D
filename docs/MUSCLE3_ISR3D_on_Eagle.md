MUSCLE3 installation on Eagle
=============================

1. Request a node with vecma account: 
-------------------------------------

```
srun --time=01:00:00 --account=vecma2020 --partition=fast --pty /bin/bash
```

Note: 
a) has to has account association, otherwise the queue priority is 0. 
b) ``partition=fast`` (1 hour operation time) is faster than ``partition=standard`` (1 week operation time).


2. Install MUSCLE3  python (anaconda recommened)
------------------------------------------------

```
pip install muscle3
```

3. Install MUSCLE3 c++ version:
-------------------------------

a) ``git clone https://github.com/multiscale/muscle3.git``

b) ``MUSCLE_ENABLE_MPI=1 make``

c) ``PREFIX=~/muscle3 make install`` (separate fold recommended)

Note: 

a) Load gcc and openmpi module before compiling. On Eagle, just use ``module load openmpi/4.0.0_gcc620``. It will load gcc 6.2 and corresponding mpi

b) NB: an error occurs when I compile MUSCLE3 on Eagle with anaconda environment on. So exit conda env before make


4. Before running ISR3D, link MUSCLE3:
--------------------------------------
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/plgrid/plgdye/muscle3/lib
```

Install ISR3D:
==============

1. Download ISR3D from gitlab and switch to branch ``muscle3-port``
-------------------------------------------------------------------

a) ``git clone https://gitlab.computationalscience.nl/pavel/ISR3D.git``

b) ``git checkout muscle3-port``

2. Compile ISR3D model
----------------------

a) Change the directory of MUSCLE3 in build.zorin.sh file

b) activate Java 8: ``module load java8/jdk1.8.0_40``

c) ANT for Java needs to be installed separately: ``wget http://apache.cs.uu.nl//ant/binaries/apache-ant-1.10.8-bin.zip``

c) unzip the folder and export $PATH:  ``export ANT_HOME=/home/plgrid/plgdye/Download/apache-ant-1.10.8 && export PATH=${PATH}:${ANT_HOME}/bin``

d) Download latest Cmake (Eagle default 3.4.3, required 3.6.3) : ``wget https://github.com/Kitware/CMake/releases/download/v3.18.0-rc4/cmake-3.18.0-rc4.tar.gz``

e) Compile and install:  ``./configure --prefix=/home/plgrid/plgdye/cmake && make && make install``

f) Export $PATH: ``export PATH=$HOME/cmake/bin:$PATH``

g) build ISR3D: ``./build.eagle.sh``

Note:

a) gcc > 6.0.0
