### GNU compiler, generic version
#
# Test environment
#   - GCC     4.4.7
#   - OpenMPI 1.10.3
#   - LAPACK  3.6.0
#
set(TARGET_SUFFIX               ".cpu")

set(ARCH                        "-march=native")
set(OPENMP_FLAGS                "-fopenmp")
set(LAPACK_FLAGS                "-llapack -lblas")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "")

set(Fortran_FLAGS_General       "-cpp")
set(C_FLAGS_General             "-std=c99 -Wall")

set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-fno-strict-aliasing -O3")
set(CMAKE_C_COMPILER            "mpicc")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-fno-strict-aliasing -O3")


########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling with GNU compiler (GCC)")
set(CMAKE_SYSTEM_PROCESSOR "generic")
