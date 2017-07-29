!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcPlateNu.F90                                !
!    CONTAINS: subroutine CalcPlateNu                     !
!                                                         ! 
!    PURPOSE: Calculate the D(vz)Dx number at the top     !
!     and bottom plates and output to a file.             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcPlateNu
      use param
      use local_arrays, only: vz
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: j,i
      real ::  nuslow, nusupp,vz_wall
      real :: del,deln
  

      nuslow = 0.d0
      nusupp = 0.d0
!      del  = 1.0/(xc(2)-xc(1))
!      deln = 1.0/(xc(nx)-xc(nxm))

      del  = 1.0/(xc(2)-xm(1))
      deln = 1.0/(xc(nx)-xm(nxm))

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,vz,del,deln) &
!$OMP   SHARED(nxm,nx) &
!$OMP   PRIVATE(i,j) &
!$OMP   REDUCTION(+:nuslow) &
!$OMP   REDUCTION(+:nusupp)

               vz_wall=0.0
      do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
         if (slipbp(j,i).eq.2) vz_wall=vz(1,j,i)
           nuslow = nuslow + (vz(1,j,i)-vz_wall)*del
         if (sliptp(j,i).eq.2) vz_wall=vz(nxm,j,i)
           nusupp = nusupp + (vz(nxm,j,i)-vz_wall)*deln
        enddo
      end do
!$OMP END PARALLEL DO

      nuslow = nuslow / (nzm*nym)
      nusupp = nusupp / (nzm*nym)

      call MpiSumRealScalar(nuslow)
      call MpiSumRealScalar(nusupp)

   if(ismaster) then
   open(97,file="shear_plate.out",status='unknown', access='sequential',position='append')
   write(97,546) time, nuslow/ReTau, nusupp/ReTau
   close(97)
   endif
546  format(F12.7,2(1x,F12.8))

      return         
      end                                                               
