cmake_minimum_required(VERSION 2.8)

# for toolchain
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/platform)

macro (set_unless_def ENV_VAL DEFAULT_VAL)
  if (NOT DEFINED ${ENV_VAL})
    set(${ENV_VAL} ${DEFAULT_VAL})
  endif ()
endmacro (set_unless_def)

macro (add_definitions_if ENV_VAL DEF_VAL)
  if (${ENV_VAL})
    add_definitions(${DEF_VAL} ${ARGN})
  endif ()
endmacro (add_definitions_if)

macro (add_definitions_unless ENV_VAL DEF_VAL)
  if (NOT ${ENV_VAL})
    add_definitions(${DEF_VAL} ${ARGN})
  endif ()
endmacro (add_definitions_unless)


### Set target
set_unless_def(BUILD_TARGET  "ms")

if (${BUILD_TARGET} STREQUAL "ms")
  add_definitions(-DARTED_MS)
  message(STATUS "Build target is Multi-Scale (MS)")
else ()
  add_definitions(-DARTED_SC)
  message(STATUS "Build target is Single-Cell (SC)")
endif ()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()


### Project settings
project(ARTED Fortran C)

message(STATUS "Target platform ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

set(TARGET_NAME         "ARTED_${BUILD_TARGET}${TARGET_SUFFIX}")
set(CMAKE_Fortran_FLAGS "${ARCH} ${OPENMP_FLAGS} ${Fortran_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS}")
set(CMAKE_C_FLAGS       "${ARCH} ${OPENMP_FLAGS} ${C_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS}")
set(EXTERNAL_LIBS       "${LAPACK_FLAGS}")

set(MODULE_LIB  modules)
set(PREP_LIB    preparation)
set(GS_LIB      GS)
set(RT_LIB      RT)
set(COMMON_LIB  common)
set(STENCIL_LIB stencil)
set(FDTD_LIB    FDTD)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

# Bugfix cmake 2.6.x: -i_dynamic flag is deperecated when using Intel Compiler
if (CMAKE_Fortran_COMPILER_ID)
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
  endif ()
endif ()


### Set optimize options
option(USE_TLOG            "Use TLOG library"               OFF)
option(USE_PAPI            "Use Performance API"            OFF)
option(DEBUG               "Enable debug"                   OFF)
option(OLD_VERSION         "Use old code"                   OFF)
option(OPT_CURRENT         "Optimize current code with RT"  ON)
option(OPT_STENCIL         "Optimize stencil code"          ON)
option(STENCIL_WITH_C      "Use C code at stencil"          ON)
option(ARRAY_PADDING       "Use array padding"              ON)
option(LOOP_BLOCKING       "Use loop blocking"              OFF)
option(DOMAIN_POWER_OF_TWO "Stencil domain is power of two" OFF)
option(EXPLICIT_VEC        "Use explicit vectorization"     ON)
option(PROFILE_THREADS     "Profile each OpenMP Threads"    OFF)

### TLOG time profiler from CCS, Univ. TSUKUBA
if (USE_TLOG)
  add_definitions(-DARTED_USE_TLOG)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} tlog)
endif ()

### Performance API from ICL, Univ. TENNESSEE
if (USE_PAPI)
  add_definitions(-DARTED_USE_PAPI)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} papi)
endif ()

add_definitions_if(DEBUG        -DARTED_DEBUG)
add_definitions_if(OPT_CURRENT  -DARTED_CURRENT_OPTIMIZED)

add_definitions_if(ENABLE_REDUCE_FOR_MANYCORE -DARTED_REDUCE_FOR_MANYCORE)

add_definitions_if(PROFILE_THREADS -DARTED_PROFILE_THREADS)

if (OPT_STENCIL)
  add_definitions(-DARTED_STENCIL_OPTIMIZED)

  add_definitions_if(ARRAY_PADDING       -DARTED_STENCIL_PADDING)
  add_definitions_if(LOOP_BLOCKING       -DARTED_STENCIL_LOOP_BLOCKING)
  add_definitions_if(DOMAIN_POWER_OF_TWO -DARTED_DOMAIN_POWER_OF_TWO)

  if (ENABLE_STENCIL_WITH_C)
    add_definitions_if(STENCIL_WITH_C -DARTED_STENCIL_WITH_C)
  elseif (STENCIL_WITH_C)
    message(WARNING ": Stencil code with C-lang: target platform is not supported.")
    set(STENCIL_WITH_C OFF)
  endif ()

  if (ENABLE_EXPLICIT_VECTORIZATION)
    add_definitions_if(EXPLICIT_VEC -DARTED_EXPLICIT_VECTORIZATION)
  elseif (EXPLICIT_VEC)
    message(WARNING ": Explicit vectorization: target platform is not supported.")
    set(EXPLICIT_VEC OFF)
  endif ()
else ()
  add_definitions(-DARTED_STENCIL_ORIGIN)
endif ()

add_definitions("${ADDITIONAL_MACRO}")


### Build target
add_subdirectory(modules)
add_subdirectory(preparation)
add_subdirectory(GS)
add_subdirectory(RT)
add_subdirectory(common)
add_subdirectory(stencil)
add_executable(${TARGET_NAME} main/${BUILD_TARGET}.f90)

add_dependencies(${TARGET_NAME} ${MODULE_LIB} ${STENCIL_LIB})
add_dependencies(${PREP_LIB}    ${MODULE_LIB})
add_dependencies(${GS_LIB}      ${MODULE_LIB})
add_dependencies(${RT_LIB}      ${MODULE_LIB})
add_dependencies(${COMMON_LIB}  ${MODULE_LIB})
add_dependencies(${STENCIL_LIB} ${MODULE_LIB})

set(LIBRARIES
    ${MODULE_LIB}
    ${PREP_LIB}
    ${GS_LIB}
    ${RT_LIB}
    ${COMMON_LIB}
    ${STENCIL_LIB}
    )

if (${BUILD_TARGET} STREQUAL "ms")
  add_subdirectory(FDTD)
  add_dependencies(${FDTD_LIB} ${MODULE_LIB})
  set(LIBRARIES ${LIBRARIES} ${FDTD_LIB})
endif ()

target_link_libraries(${TARGET_NAME} ${LIBRARIES} ${EXTERNAL_LIBS})

get_directory_property(MACRO_LOG COMPILE_DEFINITIONS)
message(STATUS "${MACRO_LOG}")
