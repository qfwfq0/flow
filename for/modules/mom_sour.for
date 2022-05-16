      module mom_sour
    
      use defs
!!#! Marco, subroutines and functions
      contains
       subroutine get_momsour(mmomsour,momsour_shape)

!       implicit none
       integer :: p
       complex (wp) , intent(inout)  :: mmomsour(0:lx,0:2)
       integer (wp) , intent(in)  :: momsour_shape

       select case (momsour_shape)
       case (0) !#! Marco, radially constant momentum source
        do p=0,2
         do ix=0,lx
          mmomsour(ix,p) = cmplx(m0_momsour(p),0.d0)
         enddo
        enddo
       
       case (1) !#! Marco, momentum source to get a Gaussian velocity in
!#! the z direction
         write(*,*) 'momsour to get a Gaussian vz',mue
         write(*,*) 'momsour intensity',m0_momsour
         write(*,*)  gauss_maxvz0,gauss_fwhm,gauss_center
        do ix=1,lx !#! remind, momsour is on the X1 mesh
         momsour(ix,2) = (- gauss_maxvz0 * mue) / (x(ix)*gauss_fwhm**2.) 
     .    * exp(-((x(ix)-gauss_center)**2.d0)/(2.d0*gauss_fwhm**2.))
     .    * (gauss_center - 2.d0*x(ix) 
     .    + (x(ix) * (x(ix)-gauss_center)**2.d0) / (gauss_fwhm**2.))
        enddo
        momsour(0,2) = (0.d0,0.d0)
!        write(*,*) 'ins_get_mom_sour'
!        write(*,*) momsour(:,2)
!        write(*,*) '**********************'
       end select

        return
       end subroutine get_momsour

       subroutine compute_mean_flow(momsour_shape)

!       implicit none
        integer (wp) , intent(in)  :: momsour_shape
        write(*,*) 'inside_compute_mean_flow', momsour_shape

        select case (momsour_shape)
        case (0) !#! Marco, radially constant momentum source 
! in z direction only
          do ix=0,lx
!#! Marco, radially constant momentum source
           vz(ix,0) = - momsour(0,2) / (4.d0 * mue) 
     .                * (x2(ix)**2.d0 - 1.d0)
          enddo
          vz(-1,0) = vz(0,0)
        
         case (1) !#! Marco, momentum source to get a Gaussian velocity in
!          write(*,*) 'ins_com_mean'
          do ix=0,lx
!#! Marco, radially constant momentum source
           vz(ix,0) =  gauss_maxvz0 
     .              * (exp(-(x(ix) - gauss_center)**2.d0 
     .              / (2.d0*gauss_fwhm**2.)))
!           write(*,*) x(ix), vz(ix,0)
          enddo
          vz(-1,0) = vz(0,0)
!          write(*,*) '##############'
         end select

        return
       end subroutine compute_mean_flow

      end module mom_sour
