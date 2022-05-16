      module mom_sour
!!#! Marco, subroutines and functions
      contains
       subroutine get_momsour(momsour)

!       implicit none
       include 'cyl1.inc.for'
       integer :: p

         write(*,*) 'm0_momsour= ',m0_momsour,lx
         do p=0,2
          do ix=0,lx
!#! Marco, radially constant momentum source
           momsour(ix,p) = cmplx(m0_momsour(p),0.d0)
          enddo
         enddo

        return
       end subroutine get_momsour

      end module mom_sour
