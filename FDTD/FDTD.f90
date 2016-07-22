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
subroutine write_result_all
  use Global_Variables
  implicit none
  integer :: iter_t,iter_i
  integer ix_m,iy_m
  character(30) wf,citer

  do iter_i=0,Nt/Nstep_write

    if(myrank == mod(iter_i,Nprocs)) then
      iter_t=iter_i*Nstep_write
      write(citer,'(I6.6)')iter_t
      wf=trim(directory)//trim(SYSname) &
        &//'_Ac_'//trim(citer)//'.out'
      open(902,file=wf)

      if(NY_m==1) then
        do ix_m=NXvacL_m,NXvacR_m
          write(902,'(10e26.16E3)') ix_m*HX_m,data_out(1:9,ix_m,1,iter_i)
        end do
      else
        do iy_m=1,NY_m
          do ix_m=NXvacL_m,NXvacR_m
            write(902,'(11e26.16E3)') ix_m*HX_m,iy_m*HY_m,data_out(1:9,ix_m,iy_m,iter_i)
          end do
        end do
      end if
      close(902)
    end if
! 1000 format(1x,2(1pe10.3,1x))
  end do
  return
end subroutine write_result_all
!===============================================================
subroutine write_result(iter)
  use Global_Variables
  implicit none
  integer i1,i2,i3,i4,i5,i6,i7,i8
  integer iter,ix_m,iy_m
  character(30) wf
  
  i1=iter/10000
  i2=mod(iter,10000)
  i3=i2/1000
  i4=mod(i2,1000)
  i5=i4/100
  i6=mod(i4,100)
  i7=i6/10
  i8=mod(i6,10)
  wf=trim(directory)//trim(SYSname) &
       &//'_Ac_'//char(48+i1)//char(48+i3)//char(48+i5)//char(48+i7)//char(48+i8)//'.out'
  open(902,file=wf)

  if(NY_m==1) then
    do ix_m=NXvacL_m,NXvacR_m
      write(902,'(10e26.16E3)') ix_m*HX_m,Ac_new_m(2,ix_m,1),Ac_new_m(3,ix_m,1) &
        &,Elec(2,ix_m,1),Elec(3,ix_m,1),j_m(2,ix_m,1),j_m(3,ix_m,1) &
        &,energy_elemag(ix_m,1),energy_elec(ix_m,1),energy_total(ix_m,1)
    end do
  else
    do iy_m=1,NY_m
      do ix_m=NXvacL_m,NXvacR_m
        write(902,'(11e26.16E3)') ix_m*HX_m,iy_m*HY_m,Ac_new_m(2,ix_m,iy_m),Ac_new_m(3,ix_m,iy_m) &
          &,Elec(2,ix_m,iy_m),Elec(3,ix_m,iy_m),j_m(2,ix_m,iy_m),j_m(3,ix_m,iy_m) &
          &,energy_elemag(ix_m,iy_m),energy_elec(ix_m,iy_m),energy_total(ix_m,iy_m) 
      end do
    end do
  end if
  close(902)
! 1000 format(1x,2(1pe10.3,1x))
  return
end subroutine write_result
!===============================================================
subroutine write_energy(iter)
  use Global_Variables
  implicit none
  integer iter,ix_m,iy_m
  character(30) wf,citer
  
  write(citer,'(I6.6)')iter
  wf=trim(directory)//trim(SYSname) &
       &//'_energy_'//trim(citer)//'.out'
  open(903,file=wf)

  if(NY_m==1) then
    do ix_m=NXvacL_m,NXvacR_m
      write(903,'(4e26.16E3)') ix_m*HX_m,energy_elemag(ix_m,1),energy_elec(ix_m,1) &
        &,energy_total(ix_m,1)
    end do
  else
    do iy_m=1,NY_m
      do ix_m=NXvacL_m,NXvacR_m
        write(903,'(5e26.16E3)') ix_m*HX_m,iy_m*HY_m,energy_elemag(ix_m,iy_m),energy_elec(ix_m,iy_m) &
          &,energy_total(ix_m,iy_m)
      end do
    end do
  end if
  close(903)
! 1000 format(1x,2(1pe10.3,1x))
  return
end subroutine write_energy
!===============================================================
subroutine write_excited_electron(iter)
  use Global_Variables
  implicit none
  integer iter,ix_m,iy_m
  character(30) wf,citer
  
  write(citer,'(I6.6)')iter
  wf=trim(directory)//trim(SYSname) &
       &//'_exc_elec_'//trim(citer)//'.out'
  open(903,file=wf)

  if(NY_m==1) then
    do ix_m=1,NX_m
      write(903,'(2e26.16E3)') ix_m*HX_m,excited_electron(ix_m,1)
    end do
  else
    do iy_m=1,NY_m
      do ix_m=1,NX_m
        write(903,'(3e26.16E3)') ix_m*HX_m,iy_m*Hy_m,excited_electron(ix_m,iy_m)
      end do
    end do
  end if
  close(903)
! 1000 format(1x,2(1pe10.3,1x))
  return
end subroutine write_excited_electron
!===============================================================
subroutine init_Ac_ms
  use Global_variables
  implicit none
  real(8) x,y
  integer ix_m,iy_m
  real(8) Xstart
  real(8) wpulse_1
  real(8) wpulse_2
! 2D parameter  
  real(8) angle,kabs,kx,ky
  real(8) length_y

  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myrank == 0)write(*,*)'ok8-1-1'

!  BC_my='isolated'
!  BC_my='periodic'
  f0_1=5.338d-9*sqrt(IWcm2_1)  ! electric field in a.u.
  omega_1=omegaeV_1/(2d0*13.6058d0)  ! frequency in a.u.
  tpulse_1=tpulsefs_1/0.02418d0 ! pulse_duration in a.u.
  Xstart=5*HX_m
  wpulse_1=2*pi/tpulse_1
  f0_2=5.338d-9*sqrt(IWcm2_2)  ! electric field in a.u.
  omega_2=omegaeV_2/(2d0*13.6058d0)  ! frequency in a.u.
  tpulse_2=tpulsefs_2/0.02418d0 ! pulse_duration in a.u.
  wpulse_2=2*pi/tpulse_2
  T1_T2=T1_T2fs/0.02418d0 ! pulse_duration in a.u.

  aY=40000d0

  
  Elec=0d0
  Bmag=0d0
  g=0d0
  j_m=0d0
  jmatter_m=0d0
  jmatter_m_l=0d0      
      
  Ac_m=0d0
  Ac_old_m=0d0
  Ac_new_m=0d0



  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  select case(FDTDdim)
  case('1D')
     if(AE_shape == 'Esin2sin') then
     do iy_m=1,NY_m
        y=iy_m*HY_m
        do ix_m=NXvacL_m,0
           x=(ix_m-1)*HX_m
           if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
              Ac_m(3,ix_m,iy_m)=Epdir_1(3)*(-f0_1/(2*omega_1)*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1+wpulse_1))*cos((omega_1+wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1-wpulse_1))*cos((omega_1-wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light))

              
              g(3,ix_m,iy_m)=-Epdir_1(3)*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
                   &*sin(omega_1*(x+Xstart+tpulse_1*c_light)/c_light)

              Ac_m(2,ix_m,iy_m)=Epdir_1(2)*(-f0_1/(2*omega_1)*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1+wpulse_1))*cos((omega_1+wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1-wpulse_1))*cos((omega_1-wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light)) 
              
              g(2,ix_m,iy_m)=-Epdir_1(2)*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
                   &*sin(omega_1*(x+Xstart+tpulse_1*c_light)/c_light)
           endif
 
          if(x > -Xstart-(tpulse_1+T1_T2)*c_light .and. x < -Xstart-(tpulse_1+T1_T2-tpulse_2)*c_light) then
      Ac_m(3,ix_m,iy_m)=Ac_m(3,ix_m,iy_m)+Epdir_2(3)*(-f0_2/(2*omega_2)*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2+wpulse_2))*cos((omega_2+wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2-wpulse_2))*cos((omega_2-wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light))

              
              g(3,ix_m,iy_m)=g(3,ix_m,iy_m)-Epdir_2(3)*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
                   &*sin(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light)

       Ac_m(2,ix_m,iy_m)=Ac_m(2,ix_m,iy_m)+Epdir_2(2)*(-f0_2/(2*omega_2)*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2+wpulse_2))*cos((omega_2+wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2-wpulse_2))*cos((omega_2-wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light)) 

             
              g(2,ix_m,iy_m)=g(2,ix_m,iy_m)-Epdir_2(2)*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
                   &*sin(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light)
           endif
        enddo
     enddo

     do iy_m=1,NY_m
       do ix_m=NXvacL_m+1,NXvacR_m-1
         Ac_new_m(2,ix_m,iy_m)=Ac_m(2,ix_m,iy_m)+dt*g(2,ix_m,iy_m)+0.5d0*(c_light*dt/HX_m)**2 &
           *(Ac_m(2,ix_m+1,iy_m)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m-1,iy_m))+0.5d0*(c_light*dt/HY_m)**2 &
           *(Ac_m(2,ix_m,iy_m+1)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m,iy_m-1))-4*pi*dt*dt*j_m(2,ix_m,iy_m) 
         
         Ac_new_m(3,ix_m,iy_m)=Ac_m(3,ix_m,iy_m)+dt*g(3,ix_m,iy_m)+0.5d0*(c_light*dt/HX_m)**2 &
           *(Ac_m(3,ix_m+1,iy_m)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m-1,iy_m))+0.5d0*(c_light*dt/HY_m)**2 &
           *(Ac_m(3,ix_m,iy_m+1)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m,iy_m-1))-4*pi*dt*dt*j_m(3,ix_m,iy_m) 
       end do
     enddo
     
   else if(AE_shape == 'Asin2cos') then
     do iy_m=1,NY_m
       y=iy_m*HY_m
       do ix_m=NXvacL_m,0
         x=(ix_m-1)*HX_m
         if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
           Ac_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
           
           Ac_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
           
         endif
         
         if(x > -Xstart-(tpulse_1+T1_T2)*c_light .and. x < -Xstart-(tpulse_1+T1_T2-tpulse_2)*c_light ) then
           
    Ac_m(3,ix_m,iy_m)=Ac_m(3,ix_m,iy_m)-Epdir_2(3)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
       &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
          
    Ac_m(2,ix_m,iy_m)=Ac_m(2,ix_m,iy_m)-Epdir_2(2)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
                  &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
         endif
         
         x=x-dt*c_light
         if(x > -Xstart-tpulse_1*c_light+dt*c_light .and. x < -Xstart+dt*c_light) then
           Ac_new_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
           
           Ac_new_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
              
         endif
         
         if(x > -Xstart-(tpulse_1+T1_T2)*c_light+dt*c_light .and. x < -Xstart-(tpulse_1+T1_T2-tpulse_2)*c_light+dt*c_light ) then
    Ac_new_m(3,ix_m,iy_m)=Ac_new_m(3,ix_m,iy_m)&
             &-Epdir_2(3)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
             &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
           
    Ac_new_m(2,ix_m,iy_m)=Ac_new_m(2,ix_m,iy_m)&
             &-Epdir_2(2)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
             &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
         endif
       enddo
     enddo
   end if

 case('2D')

   angle=45d0*pi/180d0
   kabs=omega_1/c_light
   kx=kabs*cos(angle)
   ky=kabs*sin(angle)
   length_y=2d0*pi/ky
   HY_m=length_y/dble(NY_m)
   do iy_m=1,NY_m
     y=iy_m*HY_m
     do ix_m=NXvacL_m,0
       x=(ix_m-1)*HX_m
       if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
         Ac_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
           &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y)
         
         Ac_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
           &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y)
       endif
       if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
         Ac_new_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
            &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y-omega_1*dt)            
         
         Ac_new_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
           &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y-omega_1*dt)            
       endif
     end do
   end do
 end select




  
  if(Myrank == 0)write(*,*)'ok fdtd 02'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  select case(TwoD_shape)
  case('periodic')
    Ac_new_m(:,:,0)=Ac_new_m(:,:,NY_m)
    Ac_new_m(:,:,NY_m+1)=Ac_new_m(:,:,1)
    Ac_m(:,:,0)=Ac_m(:,:,NY_m)
    Ac_m(:,:,NY_m+1)=Ac_m(:,:,1)
    g(:,:,0)=g(:,:,NY_m)
    g(:,:,NY_m+1)=g(:,:,1)
  case('isolated')
    Ac_new_m(:,:,0)=Ac_new_m(:,:,1)
    Ac_new_m(:,:,NY_m+1)=0d0
    Ac_m(:,:,0)=Ac_m(:,:,1)
    Ac_m(:,:,NY_m+1)=0d0
    g(:,:,0)=g(:,:,1)
    g(:,:,NY_m+1)=0d0
  case default
    stop 'boundary condition is not good'
  end select
  

  return
end subroutine init_Ac_ms
!===========================================================
subroutine dt_evolve_Ac
  use Global_variables
  use timelog
  implicit none
  integer ix_m,iy_m

  call timelog_begin(LOG_DT_EVOLVE_AC)

  Ac_old_m=Ac_m
  Ac_m=Ac_new_m


!  if(iter/=0)then
    do iy_m=1,NY_m
      do ix_m=NXvacL_m+1,NXvacR_m-1
        Ac_new_m(2,ix_m,iy_m)=2d0*Ac_m(2,ix_m,iy_m)-Ac_old_m(2,ix_m,iy_m)+(c_light*dt/HX_m)**2 &
          &*(Ac_m(2,ix_m+1,iy_m)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m-1,iy_m))+(c_light*dt/HY_m)**2 &
          &*(Ac_m(2,ix_m,iy_m+1)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m,iy_m-1))-4*pi*dt*dt*j_m(2,ix_m,iy_m)

        Ac_new_m(3,ix_m,iy_m)=2d0*Ac_m(3,ix_m,iy_m)-Ac_old_m(3,ix_m,iy_m)+(c_light*dt/HX_m)**2 &
          &*(Ac_m(3,ix_m+1,iy_m)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m-1,iy_m))+(c_light*dt/HY_m)**2 &
          &*(Ac_m(3,ix_m,iy_m+1)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m,iy_m-1))-4*pi*dt*dt*j_m(3,ix_m,iy_m)
      end do
    enddo
!  else
!    do iy_m=1,NY_m
!      do ix_m=NXvacL_m+1,NXvacR_m-1
!        Ac_new_m(2,ix_m,iy_m)=Ac_m(2,ix_m,iy_m)+dt*g(2,ix_m,iy_m)+0.5d0*(c*dt/HX_m)**2 &
!          *(Ac_m(2,ix_m+1,iy_m)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m-1,iy_m))+0.5d0*(c*dt/HY_m)**2 &
!          *(Ac_m(2,ix_m,iy_m+1)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m,iy_m-1))-4*pi*dt*dt*j_m(2,ix_m,iy_m) 
!
!        Ac_new_m(3,ix_m,iy_m)=Ac_m(3,ix_m,iy_m)+dt*g(3,ix_m,iy_m)+0.5d0*(c*dt/HX_m)**2 &
!          *(Ac_m(3,ix_m+1,iy_m)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m-1,iy_m))+0.5d0*(c*dt/HY_m)**2 &
!          *(Ac_m(3,ix_m,iy_m+1)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m,iy_m-1))-4*pi*dt*dt*j_m(3,ix_m,iy_m) 
!      end do
!    enddo
!  end if
      
      
  select case(TwoD_shape)
  case('periodic')
    Ac_new_m(:,:,0)=Ac_new_m(:,:,NY_m)
    Ac_new_m(:,:,NY_m+1)=Ac_new_m(:,:,1)
  case('isolated')
    Ac_new_m(:,:,0)=Ac_new_m(:,:,1)
    Ac_new_m(:,:,NY_m+1)=0d0
  end select

  call timelog_end(LOG_DT_EVOLVE_AC)
      
  return
end subroutine dt_evolve_Ac
