  implicit none
  integer NXmin,NXmax,NYmin,NYmax,NZmin,NZmax
  integer NM,j,ix,iz,iymin,iymax,iy
  integer,allocatable :: mgrid(:,:,:)
  
!  NXmin=1; NXmax=16
!  NYmin=1; NYmax=48
!  NZmin=1; NZmax=48
  NXmin=1; NXmax=4
  NYmin=1; NYmax=10
  NZmin=1; NZmax=10
 
  allocate(mgrid(NXmin:NXmax,NYmin:NYmax,NZmin:NZmax))
!  NM=23
  NM=4
  j=0
  mgrid=0
  do ix=NXmin,NXMax
  do iz=NZmin,NM
    iymin=(NYmax/2+1)-(NM-1)+(iz-1)
	iymax=(NYmax/2+1)+(NM-1)-(iz-1)
	do iy=iymin,iymax
	  j=j+1
	  mgrid(ix,iy,iz)=1
	enddo
  enddo
  do iz=NZmax-(NM-2),NZmax
    iymin=(NYmax/2+1)-(NM-2)+(NZmax-iz)
	iymax=(NYmax/2+1)+(NM-2)-(NZmax-iz)
	do iy=iymin,iymax
	  j=j+1
	  mgrid(ix,iy,iz)=1
	enddo
  enddo
  write(*,*) j
  enddo
  
  open(7,file='macrogrid_25_4.dat')
  write(7,*) j
  do iz=NZmin,NZmax
  do iy=NYmin,NYmax
  do ix=NXmin,NXmax
    if(mgrid(ix,iy,iz) == 1) then
	  write(7,*) ix,iy,iz
	endif
  enddo
  enddo
  enddo
  close(7)
  
  stop
  end
  