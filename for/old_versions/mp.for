      module mp

!#! Marco, subroutines and functions
       contains
        subroutine apply_irrot_bc(lx,mp_percr,amp_ar,m_mpr,n_mpr
     .                           ,ph_mpr,bbr,bbt,bbz)

!         implicit none
         include 'cyl1.inc.for'
!         integer , parameter :: wp = kind(1.d0) ! working precision
!         integer , intent(in) :: lx
!         integer, intent(in)   :: m_mp, n_mp
!         real(wp) , intent(in) :: amp_ar
!         integer(wp) , intent(in) :: n_mpr, m_mpr
         real(wp) , intent(in) :: mp_percr
         real(wp) , intent(in) :: ph_mpr
         complex(wp) , intent(out) , dimension(0:lx) :: bbr
         complex(wp) , intent(out) , dimension(-1:lx) :: bbt, bbz
!         complex(wp) :: i
         real(wp)  :: dumr, dumt, dumz
         real(wp)  :: sign_r, sign_t, dx
!         real(wp)  :: rr
!         real(wp)  :: pi, pitol
         real(wp) , dimension(0:lx) :: axx
         real(wp) , dimension(-1:lx) :: axx2
         integer :: ix
!         data   pi/3.14159265358979323846d0/
!         data   i/(0.d0,1.d0)/
!         data   rr/4.d0/
!         data   pitol/1.e-4/

         write(*,*) 'prova RR module',rr,wp,lx,m_mp,n_mp,i,pi,pitol
!#! Marco, define meshes
         dx = real(1.d0/lx,wp)
         do ix=0, lx
          axx(ix) = real(ix,wp) * dx
         enddo
         do ix=-1, lx
          axx2(ix) = real((ix + 1),wp) * dx - dx/2.d0
         enddo

         sign_r = signr(ph_mpr)
         sign_t = signt(ph_mpr)
!         write(*,*) 'DIAG_apply_before',ph_mpr,dtan(ph_mpr)
!         do ix=0, lx
!           dumr = dsqrt((amp_ar * 0.5d0 
!     .           * (BESSI(m_mpr+1, n_mpr/RR*axx(IX)) 
!     .           +  BESSI(m_mpr-1,n_mpr/RR*axx(IX))))**2.d0)
!
!          if ( (abs(ph_mpr) .lt. pitol) .or. 
!     .         (abs(ph_mpr-pi) .lt. pitol) .or.
!     .         (abs(ph_mpr-2.d0*pi) .lt. pitol) ) then
!             bbr(ix) =  cmplx(0.d0 * sign_r * dumr
!     .               , (- 1.d0) * sign_r * dumr)
!          else
!             bbr(ix) = cmplx(sign_r * dumr 
!     .               , (- 1.d0) / dtan(ph_mpr) * sign_r * dumr)
!     .               * ( dsqrt( dtan(ph_mpr)**2.d0   
!     .               / (1.d0 + dtan(ph_mpr)**2.d0) ))
!          endif
!         enddo
!         write(*,*) 'DIAG_apply1',sign(1,n_mpr)
!         write(*,*) 'DIAG_apply2',sign_r * dumr
!     .              ,(- 1.d0) / dtan(ph_mpr) * sign_r * dumr

         do ix=0, lx
          dumt = dsqrt((1.d0/axx2(IX) * amp_ar * m_mpr
     .           * RR / n_mpr * BESSI(m_mpr, n_mpr/RR*axx2(IX)))**2.d0 )
          dumz = dsqrt((amp_ar 
     .           * BESSI(m_mpr, n_mpr/RR*axx2(IX)))**2.d0 )
          if ((abs(ph_mpr-pi/2.d0) .lt. pitol) .or.
     .        (abs(ph_mpr-3.d0*pi/2.d0) .lt. pitol)    ) then
           bbt(ix) = cmplx( sign(1,n_mpr) * 0.d0 * sign_t * dumt
     .            ,  sign(1,n_mpr) * sign_t * dumt )
           bbz(ix) = cmplx(  sign_t * dumz * 0.d0
     .            ,  sign_t * dumz )
          else
           bbt(ix) = cmplx( sign(1,n_mpr) * sign_t * dumt
     .            , sign(1,n_mpr) * dtan(ph_mpr) * sign_t * dumt)
     .            / dsqrt(1.d0 + dtan(ph_mpr)**2.d0) 
           bbz(ix) = cmplx(  sign_t * dumz
     .            ,  dtan(ph_mpr) * sign_t * dumz)
     .            / dsqrt((1.d0 + dtan(ph_mpr)**2.d0)) 
          endif
         enddo
!#! Marco, regularity condition
         if (m_mpr .eq. 0) then 
          bbt(-1) = -bbt(0)
          bbz(-1) = +bbz(0)
         else if (m_mpr .eq. 1) then
          bbt(-1) = +bbt(0)
          bbz(-1) = -bbz(0)
         else
          bbt(-1) = -bbt(0)
          bbz(-1) = -bbz(0)
         endif

!#! Marco, inverto i due campi perché ho un segno meno rispetto a
!PIXIE3D
!!#! Marco, 09-05-2016: provo a togliere il segno meno (rimesso)
!!#! Marco, 13-07-2016: provo a togliere il segno meno
!         bbt = - bbt
!         bbz = - bbz

!#! here I should insert a divergence cleaning part, computing br from
!bt,bz, mimando BR AUS BT +BZ
!!#! Marcp, ricorda che il br deve avere una fase pi/2 maggiore di
!quella di bt,bz
!         write(*,*) 'after'
         bbr(0) = cmplx(0.d0,0.d0)
         do ix=1, lx 
          bbr(ix) = 1.d0 / axx(ix) * (
     .              i * dx *(-m_mpr * bbt(ix-1) 
     .              - n_mpr / rr * axx2(ix-1) * bbz(ix-1) ) 
     .              + bbr(ix-1) * axx(ix-1) )
         enddo
         if (m_mpr .eq. 1) then
          bbr(0) = (4.d0*bbr(1) - bbr(2)) / 3.d0
          bbt(-1) = - bbt(0) +2.d0*i*bbr(0)
         endif

!!#! Marco, check wheter the amplitude at edge is correct
!#! Marco, not necessary
!!          write(*,*) 'DIAGcrt_fct',mp_percr
!          crt_fct = mp_percr / abs(bbr(lx))
!!          write(*,*) 'DIAGcrt_fct',crt_fct
!          bbr = bbr * crt_fct
!          bbt = bbt * crt_fct
!          bbz = bbz * crt_fct
!
!         write(*,*) 'DIAG_apply2',bbr
        end subroutine apply_irrot_bc

      function signr(phmpr)
!#! Marco, building signs of the real and imaginary parts of magnetic
!field components, depending on the value of ph_mp
!#!Marco, cambiato il 14/07/2016
       implicit none
       real (8) , intent(in) :: phmpr
       real (8) :: signr, pi
       DATA   pi/3.14159265358979323846d0/
        if (phmpr .ge. 0.d0 .and. phmpr .lt. pi/2.d0) then
         signr = - 1. 
        else if (phmpr .ge. pi/2.d0 .and. phmpr .lt. pi) then
         signr = - 1. 
        else if (phmpr .ge. pi .and. phmpr .lt. 3.d0*pi/2.d0) then
         signr = + 1. 
        else if (phmpr .ge. 3.d0*pi/2.d0 .and. phmpr .lt. 2.d0*pi) then
         signr = + 1. 
        endif
      end function signr
       
      function signt(phmpr)
!#! Marco, building signs of the real and imaginary parts of magnetic
!field components, depending on the value of phmpr
       implicit none
       real (8) , intent(in) :: phmpr
       real (8) :: signt, pi
       DATA   pi/3.14159265358979323846d0/
        if (phmpr .ge. 0.d0 .and. phmpr .lt. pi/2.d0) then
         signt = + 1. !* -1.
        else if (phmpr .ge. pi/2.d0 .and. phmpr .lt. pi) then
         signt = - 1. ! * -1.
        else if (phmpr .ge. pi .and. phmpr .lt. 3.d0*pi/2.d0) then
         signt = - 1.  !* -1.
        else if (phmpr .ge. 3.d0*pi/2.d0 .and. phmpr .lt. 2.d0*pi) then
         signt = + 1.  !* -1.
        endif
      end function signt

      function get_phase(field)
!#! Marco, get the correct phase from atan
       implicit none
       complex (8) , intent(in) :: field
       real (8) :: get_phase, pi, phase
       DATA   pi/3.14159265358979323846d0/

       phase = atan( aimag(field) / real(field) )

       if (aimag(field) .ge. 0.d0) then
        if (real(field) .ge. 0d0) then
         get_phase = phase         
        else
!#! Marco, phase in questo caso è negativa
         get_phase = phase + pi
        endif
       else 
        if (real(field) .ge. 0d0) then
!#! Marco, phase in questo caso è negativa
         get_phase = 2.d0*pi + phase
        else 
         get_phase = phase + pi
        endif
       endif

      end function get_phase
!************************************************************************
!*                                                                      *
!*    Program to calculate the first kind modified Bessel function of   *
!*    integer order N, for any REAL X, using the function BESSI(N,X).   *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*                                                                      *
!*    SAMPLE RUN:                                                       *
!*                                                                      *
!*    (Calculate Bessel function for N=2, X=0.75).                      *
!*                                                                      *
!*    Bessel function of order  2 for X =  0.7500:                      *
!*                                                                      *
!*         Y =  0.73666878E-01                                          *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!*                                                                      *
!*                               F90 Release 1.1 By J-P Moreau, Paris.  *
!*                                        (www.jpmoreau.fr)             *
!*                                                                      *
!*   Version 1.1: corected value of P4 in BESSIO (P4=1.2067492 and not  *
!*                1.2067429) Aug. 2011.                                 *
!************************************************************************
c$$$PROGRAM TBESSI
c$$$
c$$$  REAL*8  BESSI, X, Y
c$$$  INTEGER N
c$$$
c$$$  N=2
c$$$  X=0.75D0
c$$$
c$$$  Y = BESSI(N,X)
c$$$
c$$$  write(*,10) N, X
c$$$  write(*,20) Y
c$$$
c$$$  stop
c$$$
c$$$ 10    format (/' Bessel Function of order ',I2,' for X=',F8.4,':')
c$$$ 20     format(/'      Y = ',E15.8/)
c$$$
c$$$END

! ----------------------------------------------------------------------
      FUNCTION BESSI(N,X)
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      REAL(8) :: X,BESSI,TOX,BIM,BI,BIP
      IF (N.EQ.0) THEN
      BESSI = BESSI0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSI = BESSI1(X)
      RETURN
      ENDIF
      IF(X.EQ.0.D0) THEN
      BESSI=0.D0
      RETURN
      ENDIF
      TOX = 2.D0/X
      BIP = 0.D0
      BI  = 1.D0
      BESSI = 0.D0
      M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
      DO 12 J = M,1,-1
      BIM = BIP+DFLOAT(J)*TOX*BI
      BIP = BI
      BI  = BIM
      IF (ABS(BI).GT.BIGNO) THEN
      BI  = BI*BIGNI
      BIP = BIP*BIGNI
      BESSI = BESSI*BIGNI
      ENDIF
      IF (J.EQ.N) BESSI = BIP
 12    CONTINUE
      BESSI = BESSI*BESSI0(X)/BI
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------
! Auxiliary Bessel functions for N=0, N=1
      FUNCTION BESSI0(X)
      REAL(8) :: X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,
     .     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,
     .     1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     .     0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     .     0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------
      FUNCTION BESSI1(X)
      REAL(8) :: X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,
     .     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     .     0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     .     -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     .     -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------

      end module mp
