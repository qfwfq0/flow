      module equilibrium
!#! Marco, subroutines and functions
      contains
       subroutine compute_equil(x,x2,bt00,bz00,vr00)

!       implicit none
       include 'cyl1.inc.for'

!       real (wp) , intent(in) :: x, x2

! daniele, 09/05/2007
! per calcolo dell'equilibrio ohmico screw pinch come in PIXIE3D
       REAL (wp) :: ohmbt(0:LX+1), ohmbz(0:LX+1)
       REAL (wp) :: ohmbt_old(0:LX+1), ohmbz_old(0:LX+1)
       REAL (wp) :: dummy(0:LX+1)

!#! Marco, new declaration
       REAL (wp) , intent(inout) :: vr00(0:LX)
       REAL (wp) , intent(inout) :: bt00(-1:LX)
       REAL (wp) , intent(inout) :: bz00(-1:LX)

       IF (IPR.EQ.-2) THEN
          write (*,*) 'Calculating ohmic equilibrium ...'

          e0_over_eta = ALPHA
          ohmbz0 = THETA0
          dr = 1./LX

             write (*,'(a,f10.3,f10.3,f10.3)')
     .              ' e0_over_eta, ohmbz0, dr: ',
     .              e0_over_eta, ohmbz0, dr

	  !Initialize iteration (sulla mesh X1!)
!          do ix=0,LX+1
          do ix=0,LX
             ohmbt(ix) = 0.5*X(ix)*e0_over_eta
             ohmbz(ix) = ohmbz0

!             write (*,'(a,i3,e10.3,e10.3)')
!     .              ' ix, ohmbt(ix), ohmbz(ix): ',
!     .              ix, ohmbt(ix), ohmbz(ix)
             
	  enddo

	  !Iteration
	  do it=1,100

!             do ix=0,LX+1
             do ix=0,LX
                ohmbt_old(ix) = ohmbt(ix)
                ohmbz_old(ix) = ohmbz(ix)
             enddo

	    !Perform integral of B_theta
             dummy(0) = 0d0
             ohmbt(0) = 0d0
!             do ix=1,LX+1
             do ix=1,LX

                bzh = 0.5*(ohmbz(ix)+ohmbz(ix-1))
                bth = 0.5*(ohmbt(ix)+ohmbt(ix-1))

                dummy(ix) = dummy(ix-1) + 
     .                 bzh**2/(bzh**2+bth**2)*x2(ix-1)*dr
     .                 /((1.+(ALET-1.)*x2(ix-1)**BEET)**GAET)
                ohmbt(ix) = e0_over_eta*dummy(ix)/x(ix)

             enddo

            !Perform integral of Bz
             dummy(0) = 0d0
!             do ix=1,LX+1
             do ix=1,LX

                bzh = 0.5*(ohmbz(ix)+ohmbz(ix-1))
                bth = 0.5*(ohmbt(ix)+ohmbt(ix-1))

                dummy(ix) = dummy(ix-1) + bth/(bzh**2+bth**2)*dr
     .                 /((1.+(ALET-1.)*x2(ix-1)**BEET)**GAET)
!                write (*,*) x2(ix-1),(1.+(ALET-1.)*x2(ix-1)**BEET)

                ohmbz(ix) = ohmbz0*exp(-e0_over_eta*dummy(ix))

             enddo

	    !Check convergence
             err = 0.
!             do ix=0,LX+1
             do ix=0,LX
                err = err
     .                 + (ohmbt(ix)-ohmbt_old(ix))**2.
     .                 + (ohmbz(ix)-ohmbz_old(ix))**2.
             enddo
             err = sqrt(0.5*err/LX)

             write (*,*) 'Equilibrium iter =',it,' ; Error=',err

             if (err < 0.1*dr**2) exit

          enddo

          write (*,*)
          write (*,'(a,i3,a,e10.2)') ' Equilibrium converged in ',it
     $           ,' iterations with error =',err

          ! Setta i campi
!          do ix=0,LX
          do ix=0,LX-1
             bz00(ix) = 0.5*(ohmbz(ix+1)+ohmbz(ix))
             bt00(ix) = 0.5*(ohmbt(ix+1)+ohmbt(ix))
          enddo

          !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
          bz_avg = 0d0
          do ix=0,LX-1
             bz_avg = bz_avg +
     .              (bz00(ix+1)*x2(ix+1) + bz00(ix)*x2(ix))*dr
          enddo

! calcolo vr00
          E0 = e0_over_eta*ETA0
          do ix=0,LX
             vr00(ix) = -E0*ohmbt(ix)/(ohmbz(ix)**2.+ohmbt(ix)**2.)
          enddo
          VRLX = vr00(LX)

          write (*,*) 
          write (*,*) 'Pinch equilibrium features:'
          write (*,*) '  Bz(r=0)=',ohmbz(0)
          write (*,*) '  Bth(r=a)=',ohmbt(LX)
          write (*,*) '  Bz avg=',bz_avg
          write (*,*) '  Theta  =',ohmbt(LX)/bz_avg
          write (*,*) '  I0     =',2.*pi*x(LX)*ohmbt(LX)
          write (*,*) '  E0/eta =',e0_over_eta
          write (*,*) '  T_flux =',pi*x(LX)**2.*bz_avg
          call flush(6)

       endif

       DO 435 IX=0,LX-1
! Daniele, 03/10/2007
! Caso Tokamak alla Coelho
! q(r) = q0 + (q1-q0)r^2
! NB qui uso C1 come q1
          IF (IPR.EQ.-3) THEN
             bz00(IX) = B0                                                
             bt00(IX) = X2(IX)*bz00(IX)/(RR*
     .              (Q0 + (C1-Q0)*X2(IX)**2.))
             vr00(IX) = 0.

! Daniele, 09/05/2007
! inserisco il caso equilibrio screw pinch ohmico
          ELSE IF (IPR.EQ.-2) THEN

! Caso Tokamak con Bz costante
! e Bt definito dal profilo di q come
! Bt = B0 / R * r / q(r)
          ELSE
           CALL EQUIL(X2(IX),bz00(IX),bt00(IX),P01(IX))                    
! Daniele, 16/05/2007
           vr00(IX) = 0.
          END IF
          BZ(IX,J0) = bz00(IX)
          BT(IX,J0) = bt00(IX)
          VR(IX,J0) = vr00(IX)
 1002     FORMAT (1X,7E12.3)                                               
 435   CONTINUE                                                          
       bz00(-1) =  bz00(0)                                              
       bt00(-1) = -bt00(0)
! Daniele, 03/10/2007
! Caso Tokamak alla Coelho
! q(r) = q0 + (q1-q0)r^2
! NB qui uso C1 come q1
       IF (IPR.EQ.-3) THEN
          BZWALL = B0                                                
          BTWALL = X(LX)*B0/(RR*
     .           (Q0 + (C1-Q0)*X(LX)**2.))
          vr00(LX) = 0.

! Daniele, 09/05/2007
! inserisco il caso equilibrio screw pinch ohmico
       ELSE IF (IPR.EQ.-2) THEN
          bzwall = ohmbz(LX)
          btwall = ohmbt(LX)
       ELSE
! fine aggiunta
          CALL EQUIL(X(LX),BZWALL,BTWALL,P0WALL)
! Daniele, 16/05/2007
          vr00(LX) = 0.
! aggiunta
       END IF
       bz00(LX) = 2.*BZWALL - bz00(LX-1)                                
       bt00(LX) = 2.*BTWALL - bt00(LX-1)                                
       P01(LX) = 2.*P0WALL - P01(LX-1)                                
       BZ(LX,J0) = bz00(LX)                                            
       BT(LX,J0) = bt00(LX)                                            
       VR(LX,J0) = vr00(LX)
       BZ(-1,J0) = BZ(0,J0)                                           
       BT(-1,J0) = -BT(0,J0)                                          
                                                                        
       A = DT * BETA * BZ(0,J0)                                       
       A22= A ** 2.                                                   

       WRITE (6,*) ' EQUILIBRIO INIZIALE FATTO: '
!       WRITE (6,*) ' BZWALL = ',BZWALL,' BTWALL = ',BTWALL
!     .        , ' POWALL = ',P0WALL
       WRITE (6,*) ' bz00(LX) = ',bz00(LX),' bt00(LX) = ',bt00(LX)
       call flush(6)
                                                                        
       end subroutine compute_equil
      end module equilibrium
