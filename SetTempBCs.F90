!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetTempBCs.F90                                 !
!    CONTAINS: subroutine SetTempBCs                      !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     temperature boundary conditions at the plates       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetTempBCs
      use param
      implicit none
      integer :: ic,jc

      do ic=1,nzm
       do jc=1,nym
        temptp(jc,ic)=0.d0
        tempbp(jc,ic)=0.d0
       enddo
      enddo
! Everywhere no-slip
      do ic=1,nzm
       do jc=1,nym
        sliptp(jc,ic)=1 ! 1 mean no-slip and 2 means slip
        slipbp(jc,ic)=1 ! 1 mean no-slip and 2 means slip// Not used as
       ! of now
       enddo
      enddo

      return
      end
!
