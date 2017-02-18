set(SOURCES
    arted_sc.f90
    arted_ms.f90
    )

add_library(${CONTROL_LIB} STATIC ${SOURCES})
