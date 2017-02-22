set(SOURCES
    dt_evolve.f90
    dt_evolve_hpsi.f90
    current.f90
    current0.f90
    k_shift_wf.f90
    hpsi_RT.f90
    init_Ac.f90
    Fourier_tr.f90
    )

add_library(${RT_LIB} STATIC ${SOURCES})
