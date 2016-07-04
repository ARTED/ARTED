###
### Compile fortran module sources
###
set(SOURCES
    global_variables_${BUILD_TARGET}.f90
    timelog.f90
    env_variables.f90
    opt_variables.f90
    env_variables_internal.c
    papi_wrap.c
    )
add_library(${MODULE_LIB} STATIC ${SOURCES})
