!
!  Copyright 2016 ARTED developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
#define TIMELOG_BEG(id) call timelog_thread_begin(id)
#define TIMELOG_END(id) call timelog_thread_end(id)

subroutine dt_evolve_hpsi
  use Global_Variables
  use timelog
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  use omp_lib
  use opt_variables
  implicit none
  integer    :: tid
  integer    :: ikb,ik,ib,i
  integer    :: iexp
  complex(8) :: zfac(4)
#ifdef ARTED_LBLK
  integer    :: idx_b,idx
  integer    :: ikb_s,ikb_e
  integer    :: ikb0,ikb1,num_ikb1
#endif
#ifdef ARTED_SC
  integer    :: loop_count
  loop_count = 0
#endif

  zfac(1)=(-zI*dt)
  do i=2,4
    zfac(i)=zfac(i-1)*(-zI*dt)/i
  end do

#ifdef ARTED_USE_NVTX
  call nvtxStartRange('dt_evolve_hpsi',2)
#endif
  call timelog_begin(LOG_HPSI)

#ifndef ARTED_LBLK

#ifdef ARTED_SC
!$omp parallel private(tid) shared(zfac) firstprivate(loop_count)
#else
!$omp parallel private(tid) shared(zfac)
!$  tid=omp_get_thread_num()

!$omp do private(ik,ib,iexp)
  do ikb=1,NKB
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    call init(ztpsi(:,4,tid),zu(:,ib,ik))
    call hpsi_omp_KB_RT(ik,ztpsi(:,4,tid),ztpsi(:,1,tid))
    call hpsi_omp_KB_RT(ik,ztpsi(:,1,tid),ztpsi(:,2,tid))
    call hpsi_omp_KB_RT(ik,ztpsi(:,2,tid),ztpsi(:,3,tid))
    call hpsi_omp_KB_RT(ik,ztpsi(:,3,tid),ztpsi(:,4,tid))
    call update(zfac,ztpsi(:,:,tid),zu(:,ib,ik))

#ifdef ARTED_CURRENT_OPTIMIZED
    call current_omp_KB_ST(ib,ik,zu(:,ib,ik))
#endif
  end do
!$omp end do

#ifdef ARTED_SC
  hpsi_called(tid) = loop_count
#endif

#else ! #ifdef ARTED_LBLK

#ifdef ARTED_SC
!$omp parallel private(tid,idx_b) shared(zfac) firstprivate(loop_count)
#else
!$omp parallel private(tid,idx_b) shared(zfac)
#endif
!$  tid=omp_get_thread_num()
  idx_b = tid * blk_nkb_hpsi

!$omp do private(num_ikb1,ikb_s,ikb_e,ikb,ik,ib,idx)
  do ikb0=1,NKB, blk_nkb_hpsi
    num_ikb1 = min(blk_nkb_hpsi, NKB-ikb0+1)
    ikb_s = ikb0
    ikb_e = ikb0 + num_ikb1-1

    call init_LBLK(ztpsi(:,:,idx_b),zu(:,:,:), ikb_s,ikb_e, 4)

    call hpsi_omp_KB_RT_LBLK(ztpsi(:,:,idx_b), ikb_s,ikb_e, 4,1)
    call hpsi_omp_KB_RT_LBLK(ztpsi(:,:,idx_b), ikb_s,ikb_e, 1,2)
    call hpsi_omp_KB_RT_LBLK(ztpsi(:,:,idx_b), ikb_s,ikb_e, 2,3)
    call hpsi_omp_KB_RT_LBLK(ztpsi(:,:,idx_b), ikb_s,ikb_e, 3,4)

    call update_LBLK(zfac,ztpsi(:,:,idx_b),zu(:,:,:), ikb_s,ikb_e)
#ifdef ARTED_CURRENT_OPTIMIZED
    call current_omp_KB_ST_LBLK(zu(:,:,:), ikb_s,ikb_e)
#endif

#ifdef ARTED_SC
    loop_count = loop_count + num_ikb1
#endif
  end do
!$omp end do

#ifdef ARTED_SC
  hpsi_called(tid) = loop_count
#endif

#endif ! ARTED_LBLK
!$omp end parallel

  call timelog_end(LOG_HPSI)
#ifdef ARTED_USE_NVTX
  call nvtxEndRange()
#endif

contains
  subroutine init(tpsi,zu)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timelog
    implicit none
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
    integer :: ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_INIT)
!dir$ vector aligned
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      tpsi(iz,iy,ix)=zu(iz,iy,ix)
    end do
    end do
    end do
    TIMELOG_END(LOG_HPSI_INIT)
  end subroutine

  subroutine update(zfac,tpsi,zu)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timelog
    implicit none
    complex(8) :: zfac(4)
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1,4)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
    integer :: ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_UPDATE)
!dir$ vector aligned
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      zu(iz,iy,ix)=zu(iz,iy,ix)+zfac(1)*tpsi(iz,iy,ix,1) &
      &                        +zfac(2)*tpsi(iz,iy,ix,2) &
      &                        +zfac(3)*tpsi(iz,iy,ix,3) &
      &                        +zfac(4)*tpsi(iz,iy,ix,4)
    end do
    end do
    end do
    TIMELOG_END(LOG_HPSI_UPDATE)
  end subroutine

#ifdef ARTED_LBLK
  subroutine init_LBLK(tpsi,zu, ikb_s,ikb_e,ia)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timelog
    implicit none
    integer :: ikb_s,ikb_e, ia
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1, 4, ikb_s:ikb_e)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1, NBoccmax, NK_s:NK_e)
    integer :: ikb,ik,ib, ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_INIT)
    do ikb=ikb_s,ikb_e
      ik=ik_table(ikb)
      ib=ib_table(ikb)
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        tpsi(iz,iy,ix, ia,ikb)=zu(iz,iy,ix, ib,ik)
      end do
      end do
      end do
    end do
    TIMELOG_END(LOG_HPSI_INIT)
  end subroutine

  subroutine update_LBLK(zfac,tpsi,zu, ikb_s,ikb_e)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timelog
    implicit none
    integer :: ikb_s,ikb_e
    complex(8) :: zfac(4)
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1, 4, ikb_s:ikb_e)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1, NBoccmax, NK_s:NK_e)
    integer :: ikb,ik,ib, ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_UPDATE)
    do ikb=ikb_s,ikb_e
      ik=ik_table(ikb)
      ib=ib_table(ikb)
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        zu(iz,iy,ix, ib,ik)=zu(iz,iy,ix, ib,ik) &
          &                +zfac(1)*tpsi(iz,iy,ix,1, ikb) &
          &                +zfac(2)*tpsi(iz,iy,ix,2, ikb) &
          &                +zfac(3)*tpsi(iz,iy,ix,3, ikb) &
          &                +zfac(4)*tpsi(iz,iy,ix,4, ikb)
      end do
      end do
      end do
    end do
    TIMELOG_END(LOG_HPSI_UPDATE)
  end subroutine
#endif
end subroutine

