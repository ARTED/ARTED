# ARTED: Ab-initio Real-Time Electron Dynamics Simulator

ARTED (Ab-initio Real-Time Electron Dynamics simulator) is an open-source
computer codes for first-principles calculations of electron dynamics and
light-matter interactions. It is based on time-dependent density functional theory
solving time-dependent Kohn-Sham equation in real time using pseudopotentials
and real-space grid representation.

ARTED has been developed in such a way that it runs
optimally in the following supercomputer platforms:

- K-computer
- Fujitsu FX100 supercomputer system
- Linux PC Cluster with x86-64 CPU
- Linux PC Cluster with Intel Xeon Phi (Knights Corner)

ARTED has been developed by ARTED developers with support from
Center for Computational Sciences, University of Tsukuba.


## Build

We use [CMake](https://cmake.org/) cross-platform build tools.

CMake detects the below configurations automatically,

- MPI Fortran/C compiler
- OpenMP flag
- LAPACK(/BLAS) libraries

### for your computer

    $ ./build
    or
    $ mkdir build_dir
    $ cd build_dir
    $ cmake .. && make

### for COMA at CCS, University of Tsukuba

    $ ./build -t sc intel avx knc
    or
    $ mkdir build_dir
    $ cd builld_dir
    $ cmake -D CMAKE_PLATFORM_TOOLCHAIN=intel-knc .. && make
    $ cmake -D CMAKE_PLATFORM_TOOLCHAIN=intel-avx .. && make

### for K-computer at RIKEN AICS

    $ ./build -t sc fujitsu k
    or
    $ mkdir build_dir
    $ cd builld_dir
    $ cmake -D CMAKE_PLATFORM_TOOLCHAIN=fujitsu-k .. && make

### for FX100 system at Nagoya University

    $ ./built fujitsu fx100
    or
    $ mkdir build_dir
    $ cd builld_dir
    $ cmake -D CMAKE_PLATFORM_TOOLCHAIN=fujitsu-fx100 .. && make

### Generic version (GNU compiler)

    $ ./built gnu generic
    or
    $ mkdir build_dir
    $ cd builld_dir
    $ cmake -D CMAKE_PLATFORM_TOOLCHAIN=gnu-generic .. && make


## Execution/Simulation

Please read to jobscript directory files.

    $ mpirun -np $NUM_MPI_PROCS ./bin/ARTED_sc.cpu < ./data/input_sc.dat


## Test Environment

### COMA at CCS, University of Tsukuba

1. Intel Compiler version 16.0.2
2. Intel MPI 5.1.3
3. Intel MKL 11.3.2

### K-computer at RIKEN AICS

1. Fujitsu Compiler version K-1.2.0-20-1

### FX100 system at Nagoya University

1. Fujitsu Compiler Driver Version 2.0.0

### GNU Compiler (Generic version)

1. GCC version 4.4.7
2. OpenMPI 1.10.3
3. LAPACK 3.6.0


## License

ARTED is available under Apache License version 2.0.

    Copyright 2016 ARTED developers
    
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
    
       http://www.apache.org/licenses/LICENSE-2.0
    
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

