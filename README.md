# ARTED: Ab-initio Real-Time Electron Dynamics Simulator

ARTED (Ab-initio Real-Time Electron Dynamics simulator) is an open-source
computer codes for first-principles calculations of electron dynamics and
light-matter interactions. It is based on time-dependent density functional theory
solving time-dependent Kohn-Sham equation in real time using pseudopotentials
and real-space grid representation.

ARTED has been developed in such a way that it runs
optimally in the following supercomputer platforms:

- K-computer
- Linux PC Cluster with x86-64 CPU
- Linux PC Cluster with Intel Xeon Phi (Knights Corner)

ARTED has been developed by ARTED developers with support from
Center for Computational Sciences, University of Tsukuba.


## Build Example

### for COMA at CCS, University of Tsukuba
    $ cd <ARTED Root Directory>
    $ module load intel intelmpi mkl
    $ ./build -t sc intel avx knc
    $ cp jobscript/template_symmetric.sh jobscript/run_symmetric.sh
    $ chmod +x jobscript/run_symmetric.sh
    $ sbatch jobscript/run_symmetric.sh
    $ ls data/

### for K-computer at RIKEN AICS
    $ cd <ARTED Root Directory>
    $ ./build -t sc fujitsu k
    $ cp jobscript/template_K.sh jobscript/run_K.sh
    $ chmod +x jobscript/run_K.sh
    $ pjsub jobscript/run_K.sh
    $ ls data/

### for GNU compiler (Generic version)
    $ cd <ARTED Root Directory>
    $ ./build -t sc gnu generic
    $ cp jobscript/template-generic.sh jobscript/run.sh
    $ chmod +x jobscript/run.sh
    $ jobscript/run.sh
    $ ls data/

## Test Environment

### COMA at CCS, University of Tsukuba

1. Intel Compiler version 16.0.2
2. Intel MPI 5.1.3
3. Intel MKL 11.3.2

### K-computer at RIKEN AICS

1. Fujitsu Compiler version K-1.2.0-19

### GNU Compiler (Generic version)

1. GCC version 4.4.7
2. OpenMPI 1.10.3
3. LAPACK 3.6.0

## License

Copyright (C) 2016  ARTED developers

Licensed under GNU GENERAL PUBLIC LICENSE (GPL) version 3.

