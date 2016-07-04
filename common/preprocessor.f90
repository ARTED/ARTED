!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of preprocessor.f90.
!
!  preprocessor.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  preprocessor.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with preprocessor.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine print_optimize_message
  implicit none
  print *, 'Preprocessor: '
#ifdef ARTED_SC
  print *, '  ARTED_SC'
#endif
#ifdef ARTED_MS
  print *, '  ARTED_MS'
#endif
#ifdef ARTED_DEBUG
  print *, '  ARTED_DEBUG'
#endif
#ifdef ARTED_USE_TLOG
  print *, '  ARTED_USE_TLOG'
#endif
#ifdef ARTED_USE_PAPI
  print *, '  ARTED_USE_PAPI'
#endif
#ifdef ARTED_CURRENT_OPTIMIZED
  print *, '  ARTED_CURRENT_OPTIMIZED'
#endif
#ifdef ARTED_STENCIL_ORIGIN
  print *, '  ARTED_STENCIL_ORIGIN'
#endif
#ifdef ARTED_STENCIL_OPTIMIZED
  print *, '  ARTED_STENCIL_OPTIMIZED'
#endif
#ifdef ARTED_STENCIL_WITH_C
  print *, '  ARTED_STENCIL_WITH_C'
#endif
#ifdef ARTED_STENCIL_PADDING
  print *, '  ARTED_STENCIL_PADDING'
#endif
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  print *, '  ARTED_STENCIL_LOOP_BLOCKING'
#endif
#ifdef ARTED_DOMAIN_POWER_OF_TWO
  print *, '  ARTED_DOMAIN_POWER_OF_TWO'
#endif
#ifdef ARTED_EXPLICIT_VECTORIZATION
  print *, '  ARTED_EXPLICIT_VECTORIZATION'
#endif
end subroutine
