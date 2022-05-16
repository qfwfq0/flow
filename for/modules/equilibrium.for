      module equilibrium

      use defs
!#! Marco, subroutines and functions
      contains
       subroutine compute_equil(ex,ex2,bt00,bz00,vr00,bzww,btww)

!       implicit none
!       include 'cyl1.inc.for'

       real (wp) , intent(in) :: ex(0:lx), ex2(0:lx)

! daniele, 09/05/2007
! per calcolo dell'equilibrio ohmico screw pinch come in PIXIE3D
       REAL (wp) :: ohmbt(0:LX+1), ohmbz(0:LX+1)
       REAL (wp) :: ohmbt_old(0:LX+1), ohmbz_old(0:LX+1)
       REAL (wp) :: dummy(0:LX+1)
       real (wp) :: e0_over_eta, ohmbz0, dr
 
!#! Marco, new declaration
       REAL (wp) , intent(inout) :: vr00(0:LX)
       REAL (wp) , intent(inout) :: bt00(-1:LX)
       REAL (wp) , intent(inout) :: bz00(-1:LX)
       real (wp) , intent(out) :: bzww,btww

       IF (IPR.EQ.-2) THEN
          write (*,*) 'Calculating ohmic equilibrium ...'

          e0_over_eta = ALPHA
          ohmbz0 = THETA0
!          dr = dble(1.d0/LX)
          dr = ex(5)-ex(4)
          write(*,*) 'dx=',dr

!             write (*,'(a,f10.3,f10.3,f10.3)')
!     .              ' e0_over_eta, ohmbz0, dr: ',
!     .              e0_over_eta, ohmbz0, dr

	  !Initialize iteration (sulla mesh X1!)
!          do ix=0,LX+1
          do ix=0,LX
             ohmbt(ix) = 0.5*X(ix)*e0_over_eta
             ohmbz(ix) = ohmbz0

!             write (*,'(a,i3,e10.3,e10.3,e10.3)')
!     .              ' ix, ohmbt(ix), ohmbz(ix): ',
!     .              ix, X(iX),ohmbt(ix), ohmbz(ix)
             
	  enddo

	  !Iteration
	  do it=1,10

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
     .                 bzh**2/(bzh**2+bth**2)*ex2(ix-1)*dr
     .                 /((1.+(ALET-1.)*ex2(ix-1)**BEET)**GAET)
                ohmbt(ix) = e0_over_eta*dummy(ix)/ex(ix)

             enddo

            !Perform integral of Bz
             dummy(0) = 0d0
!             do ix=1,LX+1
             do ix=1,LX

                bzh = 0.5*(ohmbz(ix)+ohmbz(ix-1))
                bth = 0.5*(ohmbt(ix)+ohmbt(ix-1))

                dummy(ix) = dummy(ix-1) + bth/(bzh**2+bth**2)*dr
     .                 /((1.+(ALET-1.)*ex2(ix-1)**BEET)**GAET)
!                write (*,*) ex2(ix-1),(1.+(ALET-1.)*ex2(ix-1)**BEET)

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
          bz00(-1) =  bz00(0)
          bt00(-1) = -bt00(0)

          bzww = ohmbz(LX)
          btww = ohmbt(LX)
!          write(*,*) 'inside equil', bzww,btww
          bz00(LX) = 2.*bzww - bz00(LX-1)
          bt00(LX) = 2.*btww - bt00(LX-1)
          !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
          bz_avg = 0d0
          do ix=0,LX-1
!             write(*,'(a,2f13.9)'),'bz0',ex2(ix),bz00(ix)
             bz_avg = bz_avg +
     .              (bz00(ix+1)*ex2(ix+1) + bz00(ix)*ex2(ix))*dr
          enddo

! calcolo vr00
          E0 = e0_over_eta*ETA0
          do ix=0,LX
             vr00(ix) = -E0*ohmbt(ix)/(ohmbz(ix)**2.+ohmbt(ix)**2.)
          enddo

          write (*,*) 
          write (*,*) 'Pinch equilibrium features:'
          write (*,*) '  Bz(r=0)=',ohmbz(0)
          write (*,*) '  Bth(r=a)=',ohmbt(LX)
          write (*,*) '  Bz avg=',bz_avg
          write (*,*) '  Bz,Bt(LX)=',bz00(lx),bt00(lx)
          write (*,*) '  Vr(r=a=',vr00(lx)
          write (*,*) '  Theta  =',ohmbt(LX)/bz_avg
          write (*,*) '  I0     =',2.*pi*ex(LX)*ohmbt(LX)
          write (*,*) '  E0/eta =',e0_over_eta
          write (*,*) '  T_flux =',pi*ex(LX)**2.*bz_avg
          call flush(6)

       endif

       end subroutine compute_equil

       subroutine dissipation_profiles(mmu1,mmu2,metaa,meta2,mkap)

         real (wp) , intent(inout)  :: mmu1(0:lx)
         real (wp) , intent(inout)  :: mmu2(0:lx)
         real (wp) , intent(inout)  :: metaa(0:lx)
         real (wp) , intent(inout)  :: meta2(0:lx)
         real (wp) , intent(inout)  :: mkap(0:lx)
        WRITE(6,*) ' PROFILI DEI COEFFI. DIFFUS.'                      

      ! daniele, questo è quello che c'era prima
      !            MU2(0) = MUE*(1. + (ALMU-1.)*X2(0)**BEMU)                      
      !            DO 366 IX=1,LX-1
      ! daniele, 20/06/2008
      ! inserisco parametro lambda per fare viscosità decrescente


      ! Nicholas 01/03/2022: ripristino viscosità con il coseno  
      ! Nicholas 07/07/2021: sostituisco temporaneamente il calcolo della
      ! viscosità, per aggiungere il coseno avere viscosità tendente a
      ! zero nell'edge.
      ! Nicholas 08/07/2020: conclusa la simulazione, ripristino la
      ! versione precedente.
      
        if (cos_viscos) then
              do ix=0,lx 
                MMU1(IX) = MUE*
     .             (1. + (ALMU-1.)*X(IX)**BEMU*cos(1.6*X(IX)))
                MMU2(IX) = MUE*
     .             (1. + (ALMU-1.)*X2(IX)**BEMU*cos(1.6*X2(IX)))
              enddo
        else
             do ix=0,lx
               MMU1(IX) = MUE*
     .            (1. + (ALMU-1.)*X(IX)**BEMU)**lambda
               MMU2(IX) = MUE*
     .            (1. + (ALMU-1.)*X2(IX)**BEMU)**lambda
             enddo
        endif
             do ix=0,lx
      ! daniele,17/02/2009
      ! siccome uso questo per fare anche il tokamak, avrei bisogno di gaet
                METAA(IX) = ETA0*(1. + (ALET-1.)*X(IX)**BEET)**GAET
                META2(IX) = ETA0*(1. + (ALET-1.)*X2(IX)**BEET)**GAET
!#! Marco, 31/08/2018
!#! creo strato a resistività costante, mimando uno strato di
!#! vuoto, per simulazioni tokamak
!                if (metaa(ix)/eta0 .gt. 100.) then
!                 metaa(ix) = eta0 * 100.
!                endif
!                if (meta2(ix)/eta0 .gt. 100.) then
!                 meta2(ix) = eta0 * 100.
!                endif
             enddo
       
       DO  IX=0,LX-1                                                  
          MKAP(IX) = DT*KAPPA*X2(IX)**2*(1.0-X2(IX))**2
       enddo

       end subroutine dissipation_profiles

       subroutine initial_pert
! daniele, 20/01/2006
! Inserisco flag INOPERT che è settato a 1 per saltare perturbazioni
       write (6,*) 'inopert = ',inopert
       call flush(6)
       !       IM0=1
       !       IM1=MY
       IM0=1
       IM1=1
!c se sim 2d inizializzato solo il modo m=m2d
       IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
         IM0=1
!#!Marco         IM1=5
         im1 = my
       ENDIF
       ! Daniele, 20/02/2009
       ! per perturbare solo alcuni modi
!#!Marco, scommentare per perturbazione classica dei modi
       if (.not. fasi_sganciate) then
       DO 113 IM=IM0,IM1
       write (6,*) 'Perturbo im',im,mm(im)
            M=MM(IM)
       ! Daniele, 20/02/2009
       ! per perturbare solo alcuni modi, bisogna dargli lo JANZ!
         in0=-10
         in1=-8
         IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
!#! Marco
!          IN0=N2D
!          IN1=N2D
          IN0=N2D*m
          IN1=N2D*m
         ENDIF
         DO 113 N=IN0,IN1
          IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
       	   WRITE(6,*) 'perturbo il modo',M,n2d*m
          else
       	   WRITE(6,*) 'perturbo il modo',M,N
          endif
       	  IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
!#! Marco, 19 giugno 2018
           if (tok_rot) then
            if (mod(m,m2d) .eq. 0) then
            j = janz(m,int(m/m2d*n2d))
       	    write(*,*) 'tok pert il modo',m,m/m2d*n2d
            endif
           else
            J = JANZ(M,int(n2d*m))
           endif
       	  else
       	   J = JANZ(M,N)
       	  endif
!       	  write(*,*) 'j del modo perturbato',j,lx
          J1 = JANZ(-M,-N)
          IF(J .EQ. J0) GO TO 113
! daniele,10/12/2008
! tengo solo la pert in vr, come in pixie e la metto come in pixie
          VR(0,J) = (0.,0.)
          DO 36 IX=1,LX
          call flush(6)
!#! Marco, 12/11/2014, perturbazione paradigmatica 2d
           VR(IX,J) =  cmplx(0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
     .                       0.d0)
           VR(IX,J1) = CONJG(VR(IX,J))
36           CONTINUE
113         CONTINUE
!#! end if (.not. fasi_sganciate)
       endif

!#! Marco, per perturbare vari modi con tre set di fasi diverse (a
!#!seconda del valore di n)
       if (fasi_sganciate) then
        IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
         im0 = 1
         im1 = my/2
         write(*,*) '2d computation',m2d,n2d,im0,im1
        endif
!#! Marco, tokamak case
        if (tok_rot) then
         im0 = 1
         im1 = my
        endif
        DO IM=IM0,IM1
         write (6,*) 'Perturbo im',im,mm(im)
         M=MM(IM)
       ! Daniele, 20/02/2009
       ! per perturbare solo alcuni modi, bisogna dargli lo JANZ!
         IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
!#! Marco
          if (mod(im,m2d) .eq. 0) then
           IN0=N2D*im/m2d
           IN1=N2D*im/m2d
          DO 1113 N=IN0,IN1
           if (i2d .eq. 1) then 
            write(*,*) 'perturbo 1113 caso 2d', im,n
           endif
          J = JANZ(M,N)
          call flush(6)
          J1 = JANZ(-M,-N)
          IF(J .EQ. J0) GO TO 1113
!#! Marco, fasi sganciate per simulazione 2D
            vr(0,j) = (0.d0,0.d0)
!            write(*,*) 'mod(-n,3), mod(-n,5): ',mod(-n,3),mod(-n,7)
            do ix=1, lx
             vr(ix,j) =  cmplx(mod(-n,3)*0.5*EPS1/X(IX)
     .                   *sin(pi*X(IX))**3.,
     .                   mod(-n,7)*0.5*EPS1/X(IX)
     .                   *sin(pi*X(IX))**3.)
!             if (j .eq. 1) then 
!              write(*,*) 'ix,vr(ix)',x(ix),vr(ix,j)
!             endif
            enddo 
 1113         CONTINUE
          endif
         else
         in0 = -27
         in1 = -1
 !#! Marco, tokamak case
         if (tok_rot) then
          in0 = nanf(m) 
          in1 = nanf(m) + nz(m) -1
          write(*,*) 'TOK_PERT',tok_rot,in0,in1,my
         endif
         DO 1114 N=IN0,IN1
           WRITE(6,'(a,3i3)') 'perturbo 1114 il modo',M,N,mod(-n,3)
          J = JANZ(M,N)
          call flush(6)
          J1 = JANZ(-M,-N)
          IF(J .EQ. J0) GO TO 1114
!#! Marco, fasi sganciate per simulazione 2D
                VR(0,J) = (0.,0.)
                DO 1137 IX=1,LX
                  if (mod(-n,3) .eq. 0) then
!                   if (ix .eq. 50) then 
!      write(*,*) 'n mod 3: j del modo perturbato',n,j,mod(-n,3)
!                   endif
             VR(IX,J) =  cmplx(2./3.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
     .                   1./3.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.)
                  else if (mod(-n+1,3) .eq. 0) then
!                   if (ix .eq. 50) then 
!      write(*,*) 'n+1 mod 3: j del modo perturbato',n, j,mod(-n+1,3)
!                   endif
             VR(IX,J) =  cmplx(3./7.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
     .                   4./7.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.)
                  else if (mod(-n+2,3) .eq. 0) then
!                   if (ix .eq. 50) then 
!      write(*,*) 'n+2 mod 3: j del modo perturbato',n, j,mod(-n+2,3)
!                   endif
             VR(IX,J) =  cmplx(3./11.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
     .                   8./11.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.)
                  endif
             VR(IX,J1) = CONJG(VR(IX,J))
 1137             CONTINUE
 1114         CONTINUE
         endif
      enddo !#! fine ciclo DO IM=IM0,IM1
      endif !#! fine if (fasi sganciate)


!             DO 1213 N=IN0-2,IN1-2,3
!                WRITE(6,*) 'perturbo 1213 il modo',M,N
!                J = JANZ(M,N)
!                write(*,*) 'j del modo perturbato',j
!                write(*,*) ' cmplx(3./7.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
!     .                         4./5.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.)'
!                call flush(6)
!                J1 = JANZ(-M,-N)
!                IF(J .EQ. J0) GO TO 113
!! daniele,10/12/2008
!! tengo solo la pert in vr, come in pixie
!! e la metto come in pixie
!                VR(0,J) = (0.,0.)
!                DO 1236 IX=1,LX
!                   VR(IX,J) =  cmplx(3./7.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
!     .                         4./7.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.)
!                 write(*,*) vr(ix,j)
!                VR(IX,J1) = CONJG(VR(IX,J))
!!#! cambio la fase anche di VT
! 1236             CONTINUE
! 1213         CONTINUE
!             DO 1313 N=IN0,IN1,3
!                WRITE(6,*) 'perturbo 1313 il modo',M,N
!                J = JANZ(M,N)
!                write(*,*) 'j del modo perturbato',j
!                write(*,*) 'cmplx(3./11.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
!     .         4./7.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.)'
!                call flush(6)
!                J1 = JANZ(-M,-N)
!                IF(J .EQ. J0) GO TO 113
!! daniele,10/12/2008
!! tengo solo la pert in vr, come in pixie
!! e la metto come in pixie
!                VR(0,J) = (0.,0.)
!                DO 1336 IX=1,LX
!              VR(IX,J) =  cmplx(3./11.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.,
!     .                      4./7.*0.5*EPS1/X(IX)*sin(pi*X(IX))**3.)
!                 write(*,*) vr(ix,j)
!                VR(IX,J1) = CONJG(VR(IX,J))
!!#! cambio la fase anche di VT
! 1336             CONTINUE
! 1313         CONTINUE

!       write(*,*) 'DEBUG-3',RR,m_mp,n_mp,janz(m_mp,n_mp)

       end subroutine initial_pert


      end module equilibrium
