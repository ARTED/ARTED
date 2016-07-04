###
### Compile stencil computation sources
###
set(CMAKE_Fortran_FLAGS "${ARCH} ${OPENMP_FLAGS} ${Fortran_FLAGS_General} ${Fortran_FLAGS_STENCIL}")
set(CMAKE_C_FLAGS       "${ARCH} ${OPENMP_FLAGS} ${C_FLAGS_General} ${C_FLAGS_STENCIL}")

### Soruces for current and total_energy
set(SOURCES_COMPILER_VEC
    F90/current_compiler_vec.f90
    F90/total_energy_compiler_vec.f90
    )
set(SOURCES_EXPLICIT_VEC
    C/current_explicit_vec.c
    C/total_energy_explicit_vec.c
    )

### Sources for hpsi
set(SOURCE_HPSI_ORIGINAL
    F90/hpsi_compiler_vec_original.f90
    )
set(SOURCE_HPSI_F90_COMPILER_VEC
    F90/hpsi_compiler_vec.f90
    )
set(SOURCE_HPSI_C_COMPILER_VEC
    C/hpsi_compiler_vec.c
    )
set(SOURCE_HPSI_EXPLICIT_VEC
    C/hpsi_explicit_vec.c
    )


###
if (OPT_STENCIL)
  if (STENCIL_WITH_C)
    if (EXPLICIT_VEC)
      set(SOURCES ${SOURCE_HPSI_EXPLICIT_VEC})
    else ()
      set(SOURCES ${SOURCE_HPSI_C_COMPILER_VEC})
    endif ()
  else ()
    set(SOURCES ${SOURCE_HPSI_F90_COMPILER_VEC})
  endif ()
else ()
  set(SOURCES ${SOURCE_HPSI_ORIGINAL})
endif ()

if (EXPLICIT_VEC)
  set(SOURCES ${SOURCES} ${SOURCES_EXPLICIT_VEC})
else ()
  set(SOURCES ${SOURCES} ${SOURCES_COMPILER_VEC})
endif ()

add_library(${STENCIL_LIB} STATIC ${SOURCES})
