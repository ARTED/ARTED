### Intel Compiler for Knights Corner
set(TARGET_SUFFIX               ".mic")

set(ARCH                        "-mmic")
set(SIMD_SET                    "IMCI")
set(OPENMP_FLAGS                "-qopenmp")
set(LAPACK_FLAGS                "-mkl=parallel")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "-qopt-streaming-stores always -qopt-gather-scatter-unroll=4 -qopt-ra-region-strategy=block")

set(Fortran_FLAGS_General       "-fpp -nogen-interface -std90 -warn all -diag-disable 6187,6477,6916,7025,7416")
set(C_FLAGS_General             "-Wall")

set(Fortran_FLAGS_STENCIL       "-ansi-alias -fno-alias")
set(C_FLAGS_STENCIL             "-restrict -ansi-alias -fno-alias")

set(CMAKE_Fortran_COMPILER      "mpiifort")
set(CMAKE_Fortran_FLAGS_DEBUG   "-qopt-assume-safe-padding -O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-qopt-assume-safe-padding -O3")
set(CMAKE_C_COMPILER            "mpiicc")
set(CMAKE_C_FLAGS_DEBUG         "-qopt-assume-safe-padding -O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-qopt-assume-safe-padding -O3")

set(ENABLE_STENCIL_WITH_C         1)
set(ENABLE_EXPLICIT_VECTORIZATION 1)
set(ENABLE_REDUCE_FOR_MANYCORE    1)


########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Intel Knights-Corner")
set(CMAKE_SYSTEM_PROCESSOR "knc")
