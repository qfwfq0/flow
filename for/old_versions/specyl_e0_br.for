! VERSIONE CON E0 COSTANTE

      PROGRAM SPECYL

      COMPLEX bzaini,bzafin

!      INTEGER:: inopert=0,ippcd=0,i2d=0,ivac=0

! daniele per profili tokamak alla Furth, 13/09/2006
      REAL lambda

      include 'cyl1.inc.for'

! daniele,26/11/2008
! per bc come PIXIE, caso completo
      REAL E0,VRLX,B2

! daniele, 28/12/2006
! per calcolo lineare da equilibrio ohmico
       REAL BREQ(0:LX),  BTEQ(-1:LX),  BZEQ(-1:LX)   
     .       , VREQ(0:LX),  VTEQ(-1:LX),  VZEQ(-1:LX)   

! daniele, 12/10/2007
! per tenere costanti le autofunzioni 10 (effetto toroidale)
       REAL BR10(0:LX),  BT10(-1:LX),  BZ10(-1:LX)   
     .       , VR10(0:LX),  VT10(-1:LX),  VZ10(-1:LX)   
       REAL BR20(0:LX),  BT20(-1:LX),  BZ20(-1:LX)   
     .       , VR20(0:LX),  VT20(-1:LX),  VZ20(-1:LX)   

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 21
       REAL BR21(0:LX),  BT21(-1:LX),  BZ21(-1:LX)   
     .       , VR21(0:LX),  VT21(-1:LX),  VZ21(-1:LX)   

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 32
       REAL BR32(0:LX),  BT32(-1:LX),  BZ32(-1:LX)   
     .       , VR32(0:LX),  VT32(-1:LX),  VZ32(-1:LX)   

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 43
       REAL BR43(0:LX),  BT43(-1:LX),  BZ43(-1:LX)   
     .       , VR43(0:LX),  VT43(-1:LX),  VZ43(-1:LX)   

! daniele, 09/05/2007
! per calcolo dell'equilibrio ohmico screw pinch come in PIXIE3D
       REAL ohmbt(0:LX+1), ohmbz(0:LX+1)
       REAL ohmbt_old(0:LX+1), ohmbz_old(0:LX+1)
       REAL dummy(0:LX+1)
       REAL VR0(0:LX)

      EXTERNAL FCT

      REAL HELI2D
      IF (I2D.EQ.1) THEN
         HELI2D = N2D*1./M2D
      ENDIF
      
!      write(6,*) 'prova variabile non definita: ixxx = ',inopert

      include 'cyl2.inc.for'

       write(6,*) ' '
       write(6,*) ' eta0 = ', eta0,'  mue = ',mue
       write(6,*) ' lundquist = ',1.d0/eta0,'  alet ',alet
       write(6,*) ' prandtl = ',mue/eta0,'  almu ',almu
       write(6,*) ' hart = ',1.d0/(eta0*sqrt(mue/eta0))
       write(6,*) ' itend = ',itend,' dt =',dt,'  ntp= ',itend/ib
       write(6,*) '      zeit_f (tauA) =  ',itend*dt
       write(6,*) ' correzione flusso:  itor ',itor
       write(6,*) ' ipr = ', ipr, ' lx = ',lx
       write(6,*) ' numero di modi iz ',iz
       write(6,*) ' eps1 modi ',eps1
       if (i2d.eq.1) write(6,*) ' i2d,heli2d,m2d,n2d: ',
     .        i2d,heli2d,m2d,n2d
       write(6,*) ' ivac= ',ivac
       write(6,*) ' inopert= ',inopert
       write(6,*) ' inoconv= ',inoconv
       write(6,*) ' i2d= ',i2d
       write(6,*) ' ippcd= ',ippcd
       write(6,*) ' ilinear= ',ilinear
       write(6,*) ' lambda (di norma deve essere 1) = ',lambda
       write(6,*) ' c1 (mettere >0 per regione centrale calda) = ',c1
       write(6,*) ' istatic10,af10: ',istatic10,af10
       write(6,*) ' istatic21,af21: ',istatic21,af21
       write(6,*) ' istatic32,af32: ',istatic32,af32
       write(6,*) ' istatic43,af43: ',istatic43,af43
       write(6,*) ' GAET (di norma deve essere 1) = ',GAET
       write(6,*) '   '
       call flush(6)

c    parti modificate dove compare: bza*   
c    introdotto - bzaini e bzafin (in cyl1*.for)
c               - scelta b.c.
c               - ecluso controllo sul flusso
   
       IAUS = IB
 
       IG = LX/4                                                        

       PI2 = 2.0*PI                                                     
       ITP = 0                                                          
                                                                        
       J = 0                                                            
       JANZ(0,0) = 0                                                    
       J0 = JANZ(0,0)                                                   
                                                                        
       IF(N0.EQ.0) GO TO 11                                             
       DO 25 N=NANF(0),NANF(0)+NZ(0)-1                                  
       J = J+1                                                          
       JANZ(0,N) = J                                                    
       JANZ(0,-N) = -J                                                  
   25 CONTINUE                                                          
   11  DO 10 IM=1,MY                                                    
       M = MM(IM)                                                       
       DO 10 N=NANF(IM),NANF(IM)+NZ(IM)-1                               
       J =J+1                                                           
       JANZ(M,N)=J                                                      
       JANZ(-M,-N)=-J
! daniele, 6/12/2006
! per effetto toroidale
!       if (m.eq.1 .and. n.eq.0) J10=J
!       if (m.eq.1 .and. n.eq.0) write(6,*) 'j10 = ',j10
   10 CONTINUE                                                          
       DO 61 M=1,MM(MY)                                                 
       MWERT(M) = 0                                                     
       NA(M) = 0                                                        
       NE(M) = 0                                                        
       NA1(M) = 0                                                       
   61 CONTINUE                                                          
       DO 62 M=1,MY                                                     
       MWERT(MM(M )) = 1                                                
       NA(MM(M )) = NANF(M)                                             
       NE(MM(M)) = NANF(M)+NZ(M)-1                                      
       NA1(MM(M)) = NA(MM(M))                                           
   62 CONTINUE                                                          
       MWERT(0) = 1                                                     
       NA(0) = 0                                                        
       NE(0) = 0                                                        
       NA1(0) = 0                                                       
       IF(NZ(0).EQ.0) GO TO 63                                          
       NA(0) = NANF(0)                                                  
       NE(0) = -NANF(0)                                                 
   63  DO 64 M=-1,-MM(MY),-1                                            
       MWERT(M) = MWERT(-M)                                             
       NA(M) = -NE(-M)                                                  
       NE(M) = -NA(-M)                                                  
   64  CONTINUE                                                         
       DT2 = DT*0.5                                                     
       DT22 = DT*2.0                                                    
       DT222= 2.0/DT                                                    
       ADT = ALF*DT                                                     
       DTMUE = DT*MUE                                                   
       DTKAP = DT*KAPPA                                                 
       GAM1 = GAMMA-1.0                                                 
       P0 = 47.0/(180.0*Q0**2)                                          
                                                                        
      DO 444 IM=0,MY                                                    
      M= MM(IM)                                                         
      DO 444 N=NANF(IM),NANF(IM)+NZ(IM)-1                               
                                                                        
C      WRITE(7,*) 'M',M,'    N',N,'       JANZ(M,N)=',JANZ(M,N)         
      WRITE(6,*) 'M',M,'    N',N,'       JANZ(M,N)=',JANZ(M,N)         
                                                                        
444   CONTINUE                                                          
      call flush(6)

C                                                                       
C ................ RESONANCE .....................                      
C                                                                       
!       IF (IPR.GE.0) THEN
! daniele, 9/11/2006
! modifica messa per evitare casini quando faccio simulazioni 1d (IZ=0)
! NB rimane RS = 0.4 ma non mi cambia nulla! (conta solo se è zero o non 0)
       IF (IPR.GE.0.AND.IZ.GT.0) THEN
         M = MM(1)                                                      
         N = NANF(1)
         RS = -1.                                                       
         NRAD = 100                                                     
         DRAD = 1./(NRAD - 1)                                           
         RAD = 0.0                                                      
         CALL EQUIL(RAD,BZ00,BT00,P010)                                 
         F0 = M*THETA0 + N/RR*BZ00                                      
         DO 785 IRAD=2,NRAD                                             
              RAD = (IRAD - 1)*DRAD                                     
              CALL EQUIL(RAD,BZ00,BT00,P010)                            
              F01 = M/RAD*BT00 + N/RR*BZ00                              
              DUM = F0*F01
              IF (DUM.LE.0.0) THEN                                      
                   RS = RAD - F01*DRAD/(F01 - F0)                       
                   GO TO 786                                            
              END IF                                                    
              F0 = F01                                                  
 785     CONTINUE                                                       
 786     CONTINUE                                                       
       END IF                                                           
                                                                        
C ....................MESH...............................               
                                                                        
! daniele, settembre 2006
! tolta perchè dava fastidio (nel caso tokamak)
! quanto q0>1
!       IF (IPR.EQ.-1) RS = SQRT(1.0/Q0-1.0)*C1+0.05                     
       CRS = RS/XLE                                                     
       ALC1 = ALFA1                                                     
       ALC2 = ALFA2                                                     
       CX1 = C11                                                        
       CX2 = C2                                                         
       CLX = XLE                                                        
       NX1 = LX                                                         
       TAX0 = C11*TANH(ALFA1)+C2*(TANH(ALFA2*(1.0-CRS))+TANH(ALFA2*     
     1 CRS))+1.0-C11-C2                                                 
C    STARTWERTE                                                         
       X(0) = XLA                                                       
       XX(0) = XLA                                                      
       ZZZ(0) = 0.0                                                     
       DO 534 IX=0,LX
       IF(IX.EQ.0) GO TO 407                                            
       IXXX = IX                                                        
       XST = X(IX-1)                                                   
       CALL RTNI( X(IX),FUNC,DERF,FCT,XST,1.0-8,40,IER1)                
  407  XX(IX+1) = X(IX)                                                 
       XXXX= X(IX)/XLE                                                  
       XXX1 = XXXX-CRS                                                  
       ZS = LX/XLE*(CX1*ALC1/COSH(ALC1*XXXX)**2                         
     1 +CX2*ALC2/COSH(ALC2*XXX1)**2+1.0-CX1-CX2)/TAX0                   
       ZS2(IX) = ZS*0.5                                                 
       ZSQ(IX) = ZS**2                                                  
       Z2S = -LX/XLE**2*2.0*(C2*ALFA2**2*SINH(ALFA2*XXX1)               
     1 /COSH(ALFA2*XXX1)**3+C11*ALFA1**2*TANH(ALFA1*XXXX)/COSH(ALFA1*   
     1 XXXX)**2)/TAX0                                                   
       IF(IX.EQ.0) GO TO 428                                            
       ZZZ(IX) = (Z2S+ZS/X(IX))*0.5                                     
  428  ZZ(IX) = ZZZ(IX)+ZSQ(IX)                                         
       ZZ1(IX) = ZZZ(IX)-ZSQ(IX)                                        
       Z2S2(IX) = Z2S*0.5                                               
  534 CONTINUE                                                          
!       WRITE(6,519) (X(IX) ,IX=0,LX)                                   
  519 FORMAT(1X,'X=',1P10E12.4)                                         
       DO 400 IX=0,LX-1                                                 
       X2(IX) = (X(IX)+X(IX+1))*0.5                                     
       XX2(IX+1) = X2(IX)                                               
       ZS22(IX) = (ZS2(IX)+ZS2(IX+1))                                   
       ZSQ2(IX) = (ZSQ(IX)+ZSQ(IX+1))*0.5                               
       ZZZ2(IX) = ((Z2S2(IX)+Z2S2(IX+1))+ZS22(IX)/X2(IX))*0.5           
  400 CONTINUE                                                          
! daniele, 29/11/2006
      X2(LX) = X(LX) + (X(LX)-X2(LX-1))
       ZDR = ZS
       WRITE(6,*)'ZDR: ',ZDR
C       WRITE(6,*)'X,X2,ZS2,ZSQ,ZZZ,ZZ1,Z2S2,ZS22,ZSQ2,ZZZ2',ZS,'ZS'    
C       DO 720 IX=1,10                                                  
C        WRITE(6,1003) X(IX),X2(IX),ZS2(IX),ZSQ(IX),ZZZ(IX),            
C    .                 ZZ1(IX),Z2S2(IX),ZSQ2(IX),ZZZ2(IX)               
C720    CONTINUE                                                        
C1003   FORMAT(10(E10.3,1X))                                            
                                                                        
                                                                        
       WRITE (6,*) ' MESH FATTA '                                      
                                                                        
C.................. INITIAL  PERTURBATION ..........................    
                                                                        
C      J0 = JANZ(0,0)DEFINIZIONE ANTICIPATA                             
C       DO 113 IJ=2,MPLOT                                               
C       IPM = MPLO(IJ)                                                  
C       IPN = NPLO(IJ)                                                  
C       J = JANZ(IPM,IPN)                                               
C       J1 = JANZ(-IPM,-IPN)                                            
C       M = IPM                                                         
C       DO 36 IX=0,LX                                                   
CC      VR(IX,J) = EPS1*(1.0-X(IX)**2)**2                               
CC       VT(IX,J) = I*(EPS1*(1.0-6.0*X(IX)**2+5.0*X(IX)**4))            
C       IF(M.EQ.1) THEN                                                 
C       VR(IX,J) =  EPS1*(1.0-X(IX)**2)**2                              
C       ELSE                                                            
C       VR(IX,J) =  EPS1*(1.0-X(IX)**2)**2*X(IX)**(M-1)                 
C       END IF                                                          
C       VR(IX,J1) = CONJG(VR(IX,J))                                     
C   36 CONTINUE                                                         
C       DO 114 IX=0,LX-1                                                
C       VT(IX,J)=I*(X(IX+1)*VR(IX+1,J)-X(IX)*VR(IX,J))*ZS2(IX)*2.0/M    
C       VT(IX,J1) = CONJG(VT(IX,J))                                     
C  114 CONTINUE                                                         
CC     VT(-1,J) = -VT(0,J)+2.0*I*VR(0,J)                                
C      VT(-1,J) = VT(0,J)                                               
C      VT(-1,J1) =CONJG(VT(-1,J))                                       
C       VT(LX,J) = VT(LX-1,J)*X2(LX-1)/(2.0-X2(LX-1))                   
C       VT(LX,J1) = CONJG(VT(LX,J))                                     
C  113 CONTINUE                                                         
                                                                        
C...................... INITIAL EQUILIBRIUM .....................      
                                                                        
! Daniele, 09/05/2007
! implemento il calcolo dell'equilibrio screw pinch ohmico come in PIXIE3D
! qui setto alpha, se poi voglio modificare theta mi cambia B0..
       IF (IPR.EQ.-2) THEN
          write (*,*) 'Calculating ohmic equilibrium ...'

!          e0_over_eta = 3.7
!          e0_over_eta = 3.4
          e0_over_eta = ALPHA
!          ohmbz0 = 0.706 ! per profili piatti
!          ohmbz0 = 0.885 ! per profili standard
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
!                write (*,*) x2(ix-1),(1.+(ALET-1.)*x2(ix-1)**BEET)
                
                ohmbt(ix) = e0_over_eta*dummy(ix)/x(ix)

!             write (*,'(a,i3,e10.3)')
!     .              ' ix, ohmbt(ix): ',
!     .              ix, ohmbt(ix)

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

!             write (*,'(a,i3,e10.3)')
!     .              ' ix, ohmbz(ix): ',
!     .              ix, ohmbz(ix)

             enddo

!            if (theta_ppnch /= 0d0) then
!              bz_avg = 0d0
!              do i=1,nxd
!                bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr**2
!              enddo
!              bz_avg = bz_avg/xmax**2
!              bzz0 = bth0(nxd)/theta_ppnch/bz_avg*bzz0
!            endif

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
             BZ0(ix) = 0.5*(ohmbz(ix+1)+ohmbz(ix))
             BT0(ix) = 0.5*(ohmbt(ix+1)+ohmbt(ix))
          enddo

          !Check force free diagnostic
!          do ix=0,LX-1
!             dummy(ix) =  BZ0(ix)*
!     .              (ohmbz(ix+1)-ohmbz(ix))/dr
!     $              +BT0(ix)/x2(ix)*
!     .              (ohmbt(ix+1)*x(ix+1)-ohmbt(ix)*x(ix))
!          enddo
!          write (*,*)
!          write (*,*) 'Force free error=',sqrt(sum(dummy(0:LX-1)**2)/LX)

          !Check Ohm's law
!          do ix=0,LX-1
!             dummy(ix) = e0_over_eta*BZ0(ix)
!     $              -BZ0(ix)/x2(ix)*
!     .              (ohmbt(ix+1)*x(ix+1)-ohmbt(ix)*x(ix))
!     $              +BT0(ix)*
!     .              (ohmbz(ix+1)-ohmbz(ix))/dr
!     .                 /(1.+(ALET-1.)*x2(ix-1)**BEET)**lambda
!          enddo
!          write (*,*) 'Ohms law error=',sqrt(sum(dummy(0:LX-1)**2)/LX)

          !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
          bz_avg = 0d0
          do ix=0,LX-1
             bz_avg = bz_avg +
     .              (BZ0(ix+1)*x2(ix+1) + BZ0(ix)*x2(ix))*dr
          enddo

! calcolo VR0
          E0 = e0_over_eta*ETA0
          do ix=0,LX
             vr0(ix) = -E0*ohmbt(ix)/(ohmbz(ix)**2.+ohmbt(ix)**2.)
          enddo
! Daniele, 09/01/2008
! pezza per vr0
!          vr0(LX)=0.
          VRLX = VR0(LX)

          write (*,*) 
          write (*,*) 'Pinch equilibrium features:'
          write (*,*) '  Bz(r=0)=',ohmbz(0)
          write (*,*) '  Bth(r=a)=',ohmbt(LX)
          write (*,*) '  Bz avg=',bz_avg
          write (*,*) '  Theta  =',ohmbt(LX)/bz_avg
          
          write (*,*) '  I0     =',2.*pi*x(LX)*ohmbt(LX)
          write (*,*) '  E0/eta =',e0_over_eta
          write (*,*) '  T_flux =',pi*x(LX)**2.*bz_avg

          write (*,*) 'i,r,bth0,bzz0,vrr0'
          do ix=0,LX-1
             write(*,1000) ix,x(ix),ohmbt(ix),ohmbz(ix),VR0(ix)
             write(*,1001) ix,x2(ix),BT0(ix),BZ0(ix)
          enddo
          write(*,1000) LX,x(LX),ohmbt(LX),ohmbz(LX),VR0(LX)
          write (*,*) 'i,r,q0'
          do ix=0,LX
             write(*,1000) ix,x(ix),
     .              (x(ix)*ohmbz(ix))/(RR*ohmbt(ix))
          enddo
          call flush(6)

 1000     FORMAT(I3,' ',f6.3,' ',e12.5,' ',e12.5,' ',e12.5)
 1001     FORMAT(I3,' ',f6.3,' ',e12.5,' ',e12.5)

       endif

       DO 435 IX=0,LX-1

! Daniele, 03/10/2007
! Caso Tokamak alla Coelho
! q(r) = q0 + (q1-q0)r^2
! NB qui uso C1 come q1
          IF (IPR.EQ.-3) THEN
             BZ0(IX) = B0                                                
             BT0(IX) = X2(IX)*BZ0(IX)/(RR*
     .              (Q0 + (C1-Q0)*X2(IX)**2.))
             VR0(IX) = 0.

! Daniele, 09/05/2007
! inserisco il caso equilibrio screw pinch ohmico
          ELSE IF (IPR.EQ.-2) THEN

! Caso Tokamak con Bz costante
! e Bt definito dal profilo di q come
! Bt = B0 / R * r / q(r)
          ELSE IF (IPR.EQ.-1) THEN

! Daniele, 06/11/2006
! q(r) = q_0 [1 + r^(2*lambda) [ (q_1/q_0)^lambda - 1 ]]^(1/lambda)
! q_1 = q_0 * (C1^(2*lambda) + 1)^(1/lambda)
! prima di questa data avevo pensato che fosse:
! q_1 = q_0 * (C1^2 + 1)
! ma questo vale solo nel caso lambda=1!

             BZ0(IX) = B0                                                
!             BT0(IX) = X2(IX)*BZ0(IX)/((1.0+(X2(IX)/C1)**2)*RR*Q0)
! Daniele, da quando ho iniziato a usare lambda (settembre 2006)
             BT0(IX) = X2(IX)*BZ0(IX)/(RR*Q0*
     .              (1.0 + (X2(IX)/C1)**(2.*lambda) )**(1./lambda))

c          P(IX,J0) = C1**2*0.5/((1.0+(X2(IX)/C1)**2)**2                
c     1               *(RR*Q0)**2)                                      

! Daniele, 16/05/2007
             VR0(IX) = 0.
          ELSE
             CALL EQUIL(X2(IX),BZ0(IX),BT0(IX),P01(IX))                    
! Daniele, 16/05/2007
             VR0(IX) = 0.
          END IF
          BZ(IX,J0) = BZ0(IX)
          BT(IX,J0) = BT0(IX)
          VR(IX,J0) = VR0(IX)
C       WRITE (6,1002) X2(IX), BZ(IX,J0), BT(IX,J0)                     
 1002     FORMAT (1X,7E12.3)                                               
 435   CONTINUE                                                          
       BZ0(-1) =  BZ0(0)                                              
       BT0(-1) = -BT0(0)
          
! Daniele, 03/10/2007
! Caso Tokamak alla Coelho
! q(r) = q0 + (q1-q0)r^2
! NB qui uso C1 come q1
       IF (IPR.EQ.-3) THEN
          BZWALL = B0                                                
          BTWALL = X(LX)*B0/(RR*
     .           (Q0 + (C1-Q0)*X(LX)**2.))
          VR0(LX) = 0.

! Daniele, 09/05/2007
! inserisco il caso equilibrio screw pinch ohmico
       ELSE IF (IPR.EQ.-2) THEN
          bzwall = ohmbz(LX)
          btwall = ohmbt(LX)
       ELSE IF (IPR.EQ.-1) THEN
! daniele, 03/04/06
          BZWALL = B0
!            BTWALL = X(LX)*B0/((1.0+(X(LX)/C1)**2)*RR*Q0)
          BTWALL = X(LX)*B0/(RR*Q0*
     .           (1.0 + (X(LX)/C1)**(2.*lambda) )**(1./lambda))

c          P0WALL = C1**2*0.5/((1.0+(X(LX)/C1)**2)**2                
c     1               *(RR*Q0)**2)                                      

! Daniele, 16/05/2007
          VR0(LX) = 0.
       ELSE
! fine aggiunta
          CALL EQUIL(X(LX),BZWALL,BTWALL,P0WALL)
! Daniele, 16/05/2007
          VR0(LX) = 0.
! aggiunta
       END IF

c        WRITE (6,1002) X(LX), BZWALL, BTWALL, P0WALL                   
       BZ0(LX) = 2.*BZWALL - BZ0(LX-1)                                
       BT0(LX) = 2.*BTWALL - BT0(LX-1)                                
       P01(LX) = 2.*P0WALL - P01(LX-1)                                
       BZ(LX,J0) = BZ0(LX)                                            
       BT(LX,J0) = BT0(LX)                                            
! Daniele, 16/05/2007
       VR(LX,J0) = VR0(LX)
C         WRITE(6,1002) X(LX),BZ(LX,J0),BT(LX,J0)                       
       BZ(-1,J0) = BZ(0,J0)                                           
       BT(-1,J0) = -BT(0,J0)                                          
                                                                        
       A = DT * BETA * BZ(0,J0)                                       
       A22= A ** 2.                                                   

       WRITE (6,*) ' EQUILIBRIO INIZIALE FATTO: '
       WRITE (6,*) ' BZWALL = ',BZWALL,' BTWALL = ',BTWALL
     .        , ' POWALL = ',P0WALL
       WRITE (6,*) ' BZ0(LX) = ',BZ0(LX),' BT0(LX) = ',BT0(LX)
       call flush(6)
                                                                        
C...................... INITIAL  DENSITY     .....................      
                                                                        
C      DO 706 IM=0,MY                                                   
C      M = MM(IM)                                                       
C       DO 706 N=NA1(M),NE(M)                                           
C       J= JANZ(M,N)                                                    
C       J1=JANZ(-MM(IM),-N)                                             
C         DO 706 IX=0,LX                                                
C                                                                       
C         DE(IX,J) = (0.,0.)                                            
C         PR(IX,J) = (0.,0.)                                            
C         IF( J .EQ. 0 ) DE(IX,J)= (1.,0.)                              
C         DE(IX,J1) = CONJG( DE(IX,J) )                                 
C                                                                       
C         IF( J .EQ. 0) THEN                                            
C         PRRE= 2./3.                                                   
C         IF( IX .EQ. 0 ) PRR = 1.                                      
C         IF(IX.GT.0) PRR = (1.+(ALET-1.)*X(IX)**BEET)**PRR             
C         PRR = 1./PRR                                                  
C         PR(IX,J)= (1.,0.)                                             
C         PR(IX,J)= PRR * PR(IX,J) * 0.5                                
C         WRITE(6,*)'  ',IX,'   ', PRR                                  
C         ENDIF                                                         
C                                                                       
C         PR(IX,J1) = CONJG( PR(IX,J) )                                 
                                                                        
C         IF( IX .GE. 1)THEN                                            
C         AG = 3.14 * X(IX)                                             
C         AG2= 3.14 * X2(IX)                                            
C                                                                       
C         FIX  = 0.1 *(SIN(AG))**2.                                     
C         DFIX = 0.1 * PI * SIN(2.*AG2)                                 
C                                                                       
C         FIX  = 1.0                                                    
C         DFIX = 1.0 * 1.                                               
C                                                                       
C         FIX  = 0.1 * SIN(AG)                                          
C         DFIX = 0.1 * PI * COS( AG2)                                   
C         VRP(IX,J) = 0.                                                
C         VTP(IX,J) = 0.                                                
C         IF(M .EQ. 1 .AND. N .EQ. 0)THEN                               
C         VRP(IX,J) = FIX/(X(IX)*2.)                                    
C         VRP(IX,J) = FIX/ 2.                                           
C         VTP(IX,J) = I * DFIX / 2.                                     
C         IF( IX.EQ.0) WRITE(6,*) VTP(IX,J),J,'QUI'                     
C         ENDIF                                                         
C         ENDIF                                                         
C         VRP(IX,J1)= CONJG( VRP(IX,J) )                                
C         VTP(IX,J1)= CONJG( VTP(IX,J) )                                
C                                                                       
C706    CONTINUE                                                        
C...........................................................            
                                                                        
C - PROFILES OF THE DIFFUSION COEFFICIENTS                              
                                                                        
C        WRITE(6,*) ' PROFILI DEI COEFFI. DIFFUS.'                      

! Daniele, 03/10/2007
! Caso Tokamak alla Coelho
! q(r) = q0 + (q1-q0)r^2
! NB qui uso C1 come q1
          IF (IPR.EQ.-3) THEN
            DO IX=0,LX
               ETAA(IX) = ETA0*(1. +
     .                2.*(C1-Q0)*X(IX)**2./Q0**2.+
     .                (C1-Q0)**2.*X(IX)**4./Q0**2.)
               ETA2(IX) = ETA0*(1. +
     .                2.*(C1-Q0)*X2(IX)**2./Q0**2.+
     .                (C1-Q0)**2.*X2(IX)**4./Q0**2.)
               MU1(IX) = MUE*(1. +
     .                2.*(C1-Q0)*X(IX)**2./Q0**2.+
     .                (C1-Q0)**2.*X(IX)**4./Q0**2.)
               MU2(IX) = MUE*(1. +
     .                2.*(C1-Q0)*X2(IX)**2./Q0**2.+
     .                (C1-Q0)**2.*X2(IX)**4./Q0**2.)
            ENDDO
! daniele, per tokamak, 13/09/2006
         ELSE IF (IPR.EQ.-1) THEN
            DO IX=0,LX
               ETAA(IX) = ETA0*(1. + (X(IX)/C1)**(2.*lambda) )
     .                ** (1. + 1./lambda)
               ETA2(IX) = ETA0*(1. + (X2(IX)/C1)**(2.*lambda) )
     .                ** (1. + 1./lambda)
               MU1(IX) = MUE*(1. + (X(IX)/C1)**(2.*lambda) )
     .                ** (1. + 1./lambda)
               MU2(IX) = MUE*(1. + (X2(IX)/C1)**(2.*lambda) )
     .                ** (1. + 1./lambda)
            ENDDO
         ELSE
! daniele, 07/04/2008
            if (ivac.eq.4) then
               DO IX=0,LX
                  MU1(IX) = DISVAC4(X(IX),MUE,ALMU,BEMU,C1)
                  MU2(IX) = DISVAC4(X2(IX),MUE,ALMU,BEMU,C1)
               enddo
               DO IX=0,LX
                  ETAA(IX) = DISVAC4(X(IX),ETA0,ALET,BEET,C1)
                  ETA2(IX) = DISVAC4(X2(IX),ETA0,ALET,BEET,C1)
               enddo
! daniele, 16/05/2007
            else if (ivac.eq.3) then
               DO IX=0,LX
                  MU1(IX) = DISVAC3(X(IX),MUE,ALMU,BEMU,GAMU,C1)
                  MU2(IX) = DISVAC3(X2(IX),MUE,ALMU,BEMU,GAMU,C1)
               enddo
               DO IX=0,LX
                  ETAA(IX) = DISVAC3(X(IX),ETA0,ALET,BEET,GAET,C1)
                  ETA2(IX) = DISVAC3(X2(IX),ETA0,ALET,BEET,GAET,C1)
               enddo
            else if (ivac.eq.2) then
               DO IX=0,LX
                  MU1(IX) = DISVAC2(X(IX),MUE,ALMU,C1)
                  MU2(IX) = DISVAC2(X2(IX),MUE,ALMU,C1)
               enddo
               DO IX=0,LX
                  ETAA(IX) = DISVAC2(X(IX),ETA0,ALET,C1)
                  ETA2(IX) = DISVAC2(X2(IX),ETA0,ALET,C1)
               enddo
            else if (ivac.eq.1) then
               DO IX=0,LX
                  MU1(IX) = MUE*(1. + ALMU*(sin(pi*X(ix)/2.))**BEMU)
                  MU2(IX) = MUE*(1. + ALMU*(sin(pi*X2(ix)/2.))**BEMU)
               enddo
               DO IX=0,LX
                  ETAA(IX) = ETA0*(1.+ ALET*(sin(pi*X(ix)/2.))**BEET)
                  ETA2(IX) = ETA0*(1.+ ALET*(sin(pi*X2(ix)/2.))**BEET)
               enddo
            else
! daniele, questo è quello che c'era prima
!            MU2(0) = MUE*(1. + (ALMU-1.)*X2(0)**BEMU)                      
!            DO 366 IX=1,LX-1
! daniele, 20/06/2008
! inserisco parametro lambda per fare viscosità decrescente
               DO 366 IX=0,LX
                  MU1(IX) = MUE*
     .                   (1. + (ALMU-1.)*X(IX)**BEMU)**lambda
                  MU2(IX) = MUE*
     .                   (1. + (ALMU-1.)*X2(IX)**BEMU)**lambda
 366           CONTINUE
               DO 333 IX=0,LX
! daniele,17/02/2009
! siccome uso questo per fare anche il tokamak, avrei bisogno di gaet
                  ETAA(IX) = ETA0*(1. + (ALET-1.)*X(IX)**BEET)**GAET
                  ETA2(IX) = ETA0*(1. + (ALET-1.)*X2(IX)**BEET)**GAET
!                  ETAA(IX) = ETA0*(1. + (ALET-1.)*X(IX)**BEET)
!                  ETA2(IX) = ETA0*(1. + (ALET-1.)*X2(IX)**BEET)
 333           CONTINUE
! daniele, 25/06/2008
! per avere regione calda in centro
! ricordarsi di inserirla anche in scel.for
!               if (C1 .ne. 0.) then
!                  DO IX=0,LX*C1
!                     MU1(IX) = MUE*(0.5*
!     .                      (1. + (X(IX)/C1)**BEMU))**lambda
!                     ETAA(IX) = ETA0*(0.5*
!     .                      (1. + (X(IX)/C1)**BEET))
!                  enddo
!                  DO IX=0,LX*C1-1
!                     MU2(IX) = MUE*(0.5*
!     .                      (1. + (X2(IX)/C1)**BEMU))**lambda
!                     ETA2(IX) = ETA0*(0.5*
!     .                      (1. + (X2(IX)/C1)**BEET))
!                  enddo
!               endif
            endif
         END IF

         do ix=0,lx
           write(6,*) ix, x(ix), etaa(ix), mu1(ix)
     .             , (etaa(ix)*mu1(ix))**-0.5
         end do
         do ix=0,lx
           write(6,*) ix, x2(ix), eta2(ix), mu2(ix)
     .             , (eta2(ix)*mu2(ix))**-0.5
         end do
        call flush(6)

       DO 78 IX=0,LX-1                                                  
           KAP(IX) = DT*KAPPA*X2(IX)**2*(1.0-X2(IX))**2                 
   78 CONTINUE                                                          
                                                                        
       ZEIT = 0.0                                                       
       IBIB = 0                                                         
                                                                        
       IF(.NOT.IBAND) GO TO 88                                          
                                                                        
C..................... READ FROM DISK OR FRONT-END ...............      
                                                                        
       WRITE(6,*)'INIZIO LETTURA'
        call flush(6)
                                                                        
       REWIND IBEIN                                                     
        IF(IBAND1.EQ.0)READ(IBEIN) ZEIT,LY,IZZ                          
       WRITE(6,*)'ZEIT=',ZEIT,'   LY',LY,'    IZZ',IZZ                  
                                                                        
        IF( IZZ .NE. IZ )                                               
     .  WRITE(6,*) 'NUMERO DIVERSO DI MODI ELABORATI:IZZ/IZ'

        call flush(6)
                                                                        
       DO 250 J=0,IZZ                                                   
!       DO 250 J=0,IZ
          IF(IBAND1.EQ.0)READ(IBEIN) AVT, AVZ, ABT, ABZ

! daniele, 28/12/2006
! per calcolo lineare da equilibrio ohmico
          IF (ILINEAR.EQ.1.AND.J.EQ.0) THEN
             vtEQ(-1) = AVT
             vzEQ(-1) = AVZ
             btEQ(-1) = ABT
             bzEQ(-1) = ABZ
          ENDIF

! daniele, 12/10/2007
! per tenere costanti le autofunzioni 10
          J10 = JANZ(1,0)
          J20 = JANZ(2,0)
          IF (ISTATIC10.EQ.1.AND.J.EQ.J10) THEN
             vt10(-1) = AF10*AVT
             vz10(-1) = AF10*AVZ
             bt10(-1) = AF10*ABT
             bz10(-1) = AF10*ABZ
          ENDIF
          IF (ISTATIC10.EQ.1.AND.J.EQ.J20) THEN
             vt20(-1) = AF10*AVT
             vz20(-1) = AF10*AVZ
             bt20(-1) = AF10*ABT
             bz20(-1) = AF10*ABZ
          ENDIF

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 21
          J21 = JANZ(2,-1)
          IF (ISTATIC21.EQ.1.AND.J.EQ.J21) THEN
             vt21(-1) = AF21*AVT
             vz21(-1) = AF21*AVZ
             bt21(-1) = AF21*ABT
             bz21(-1) = AF21*ABZ
          ENDIF

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 32
          J32 = JANZ(3,-2)
          IF (ISTATIC32.EQ.1.AND.J.EQ.J32) THEN
             vt32(-1) = AF32*AVT
             vz32(-1) = AF32*AVZ
             bt32(-1) = AF32*ABT
             bz32(-1) = AF32*ABZ
          ENDIF

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 43
          J43 = JANZ(4,-3)
          IF (ISTATIC43.EQ.1.AND.J.EQ.J43) THEN
             vt43(-1) = AF43*AVT
             vz43(-1) = AF43*AVZ
             bt43(-1) = AF43*ABT
             bz43(-1) = AF43*ABZ
          ENDIF

          VT(-1,J) = AVT*I                                              
          VZ(-1,J) = AVZ*I                                              
          BT(-1,J) = ABT                                                
          BZ(-1,J) = ABZ                                                
          VT(-1,-J) = CONJG(VT(-1,J))                                   
          VZ(-1,-J) = CONJG(VZ(-1,J))                                   
          BT(-1,-J) = CONJG(BT(-1,J))                                   
          BZ(-1,-J) = CONJG(BZ(-1,J))                                   
       DO 250 IX=0,LY                                                   
          IF(IBAND1.EQ.0)READ(IBEIN)     AVR, AVT, AVZ, ABR, ABT, ABZ   

! daniele, 28/12/2006
! per calcolo lineare da equilibrio ohmico
          IF (ILINEAR.EQ.1.AND.J.EQ.0) THEN
             vrEQ(ix) = AVR
             vtEQ(ix) = AVT
             vzEQ(ix) = AVZ
             brEQ(ix) = ABR
             btEQ(ix) = ABT
             bzEQ(ix) = ABZ
          ENDIF

! daniele, 12/10/2007
! per tenere costanti le autofunzioni 10
          IF (ISTATIC10.EQ.1.AND.J.EQ.J10) THEN
             vr10(ix) = AF10*AVR
             vt10(ix) = AF10*AVT
             vz10(ix) = AF10*AVZ
             br10(ix) = AF10*ABR
             bt10(ix) = AF10*ABT
             bz10(ix) = AF10*ABZ
          ENDIF
          IF (ISTATIC10.EQ.1.AND.J.EQ.J20) THEN
             vr20(ix) = AF10*AVR
             vt20(ix) = AF10*AVT
             vz20(ix) = AF10*AVZ
             br20(ix) = AF10*ABR
             bt20(ix) = AF10*ABT
             bz20(ix) = AF10*ABZ
          ENDIF

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 21
          IF (ISTATIC21.EQ.1.AND.J.EQ.J21) THEN
             vr21(ix) = AF21*AVR
             vt21(ix) = AF21*AVT
             vz21(ix) = AF21*AVZ
             br21(ix) = AF21*ABR
             bt21(ix) = AF21*ABT
             bz21(ix) = AF21*ABZ
          ENDIF

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 32
          IF (ISTATIC32.EQ.1.AND.J.EQ.J32) THEN
             vr32(ix) = AF32*AVR
             vt32(ix) = AF32*AVT
             vz32(ix) = AF32*AVZ
             br32(ix) = AF32*ABR
             bt32(ix) = AF32*ABT
             bz32(ix) = AF32*ABZ
          ENDIF

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 43
          IF (ISTATIC43.EQ.1.AND.J.EQ.J43) THEN
             vr43(ix) = AF43*AVR
             vt43(ix) = AF43*AVT
             vz43(ix) = AF43*AVZ
             br43(ix) = AF43*ABR
             bt43(ix) = AF43*ABT
             bz43(ix) = AF43*ABZ
          ENDIF

          VR(IX,J) = AVR                                                
          VT(IX,J) = AVT*I                                              
          VZ(IX,J) = AVZ*I                                              
          BR(IX,J) = ABR*I                                              
          BT(IX,J) = ABT                                                
          BZ(IX,J) = ABZ                                                
          VR(IX,-J) = CONJG(VR(IX,J))                                   
          VT(IX,-J) = CONJG(VT(IX,J))                                   
          VZ(IX,-J) = CONJG(VZ(IX,J))                                   
          BR(IX,-J) = CONJG(BR(IX,J))                                   
          BT(IX,-J) = CONJG(BT(IX,J))                                   
          BZ(IX,-J) = CONJG(BZ(IX,J))                                   

!          if (ix.eq.50) write(6,*) 'j: ',j, ' log br(50): ',
!     .           j,alog10(abs(abr))


 250   CONTINUE                                                         

! daniele, 12/10/2007
! per tenere costanti le autofunzioni 10
      IF (ISTATIC10.EQ.1.AND.AF10.EQ.0.) THEN
          DO IX=0,lX
             bz10(ix) = B0*x(ix)/(2.*RR)
!             write(6,*) 'setto ix, bz10(ix) = ',ix,bz10(ix)
          ENDDO
       ENDIF
                                                                        
C   THE PERTURBATIONS (J=-1,1) ARE NORMALIZED TO EPS1                   
       IF (INORM.EQ.1) THEN                                             
         VRMAX = 0.0                                                    
         VTMAX = 0.0                                                    
         VZMAX = 0.0                                                    
         BRMAX = 0.0                                                    
         BTMAX = 0.0                                                    
         BZMAX = 0.0                                                    
         DO 782 IX=0,LX                                                 
              YVR = ABS(REAL(VR(IX,1)))                                 
              YVT = ABS(AIMAG(VT(IX,1)))                                
              YVZ = ABS(AIMAG(VZ(IX,1)))                                
              YBR = ABS(AIMAG(BR(IX,1)))                                
              YBT = ABS(REAL(BT(IX,1)))                                 
              YBZ = ABS(REAL(BZ(IX,1)))                                 
              VRMAX = AMAX1(VRMAX,YVR)                                  
              VTMAX = AMAX1(VTMAX,YVT)                                  
              VZMAX = AMAX1(VZMAX,YVZ)                                  
              BRMAX = AMAX1(BRMAX,YBR)                                  
              BTMAX = AMAX1(BTMAX,YBT)                                  
              BZMAX = AMAX1(BZMAX,YBZ)                                  
 782     CONTINUE                                                       
         VALMAX = VRMAX                                                 
         VALMAX = AMAX1(VALMAX,VRMAX)                                   
         VALMAX = AMAX1(VALMAX,VTMAX)                                   
         VALMAX = AMAX1(VALMAX,VZMAX)                                   
         VALMAX = AMAX1(VALMAX,BRMAX)                                   
         VALMAX = AMAX1(VALMAX,BTMAX)                                   
         VALMAX = AMAX1(VALMAX,BZMAX)                                   
         FACNOR = EPS1/VALMAX                                           
         WRITE(6,1005) FACNOR                                           
 1005    FORMAT(1X,'NORMALIZING FACTOR=',1PE12.4)                       
         DO 780 J=-1,1,2                                                
         DO 781 IX=0,LX                                                 
              VR(IX,J) = VR(IX,J)*FACNOR                                
              BR(IX,J) = BR(IX,J)*FACNOR                                
 781     CONTINUE                                                       
         DO 779 IX=-1,LX                                                
              VT(IX,J) = VT(IX,J)*FACNOR                                
              VZ(IX,J) = VZ(IX,J)*FACNOR                                
              BT(IX,J) = BT(IX,J)*FACNOR                                
              BZ(IX,J) = BZ(IX,J)*FACNOR                                
 779     CONTINUE                                                       
 780     CONTINUE                                                       
! daniele, 11/05/2007
! prima era subito dopo a 780, ora lo sposto perchè avvenga sempre
! modifico anche VR0 all'equilibrio
! annullato spostamento che era dopo 88
         DO 778 IX=0,LX                                                 
            IF (IPR.EQ.-2) THEN
               VR(IX,J0) = VR0(IX)
            ELSE
               VR(IX,J0) = (0.0,0.0)
            END IF
              BR(IX,J0) = (0.0,0.0)                                     
 778     CONTINUE                                                      
         DO 776 IX=-1,LX                                                
              VT(IX,J0) = (0.0,0.0)                                     
              VZ(IX,J0) = (0.0,0.0)                                     
              BT(IX,J0) = BT0(IX)                                       
              BZ(IX,J0) = BZ0(IX)                                       
776      CONTINUE                                                       
       END IF                                                           
                                                                        
       IF (IZEIT.EQ.1) THEN                                             
         ZEIT = 0.0                                                     
         ITP = 0                                                        
       END IF                                                           
                                                                        
88     CONTINUE                                                         

       TORF1 = 0.0
       DO 991 IX=0,LX-1
          TORF1 = TORF1 + 2. * X2(IX)*REAL(BZ(IX,J0))/ZDR
991    CONTINUE
       YFITP = ( REAL(BZ(LX-1,J0))+REAL(BZ(LX,J0)) )/(2*TORF1)       
       YTH = ( REAL(BT(LX-1,J0))+REAL(BT(LX,J0)) )/(2*TORF1)       
       YQITP = X2(0)/RR*REAL(BZ(0,J0))/REAL(BT(0,J0))                 

       write(6,*) ' zeit= ',zeit,'  stato iniziale :'
       WRITE(6,*) 'x= ',X2(LX-4),' bz ',BZ(LX-4,0),' bt ',BT(LX-4,0)
       WRITE(6,*) 'x= ',X2(LX-3),' bz ',BZ(LX-3,0),' bt ',BT(LX-3,0)
       WRITE(6,*) 'x= ',X2(LX-2),' bz ',BZ(LX-2,0),' bt ',BT(LX-2,0)
       WRITE(6,*) 'x= ',X2(LX-1),' bz ',BZ(LX-1,0),' bt ',BT(LX-1,0)
       WRITE(6,*) 'x= ',X2(LX),' bz ',BZ(LX,0),' bt ',BT(LX,0)

       WRITE(6,*) '   REAL(BZ(LX-1,J0)) =',REAL(BZ(LX-1,J0)),
     .        '   REAL(BZ(LX,J0)) = ', REAL(BZ(LX,J0))
       WRITE(6,*) '   REAL(BT(LX-1,J0)) =',REAL(BT(LX-1,J0)),
     .        '   REAL(BT(LX,J0)) = ', REAL(BT(LX,J0))
       WRITE(6,*) '   0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0))) = ',
     .        0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0)))
       WRITE(6,*) '   TH =',YTH, '   F = ', YFITP
       WRITE(6,*) '   flusso toroidale= ', TORF1
       WRITE(6,*) '   Qo= ', YQITP

       call flush(6)
                                                                        
C      IF (.NOT. IBDEN) GO TO 89                                        
C      REWIND IDEIN                                                     
C      REWIND IPRIN                                                     
C      READ (IDEIN) ZEIT,LY,IZZ                                         
C      READ (IPRIN) ZEIT,LY,IZZ                                         
C      DO 710 J=0,IZZ                                                   
C       DO 710 IX=0,LY                                                  
C        READ(IDEIN) ADE                                                
C        DE(IX,J) = ADE                                                 
C        DE(IX,-J)= CONJG(DE(IX,J))                                     
C        READ(IPRIN) APR                                                
C        PR(IX,J) = APR                                                 
C        PR(IX,-J)= CONJG(PR(IX,J))                                     
C710    CONTINUE                                                        
C89     CONTINUE                                                        
                                                                        
c          t0 = second()
c          tt = t0

C                    perturbazione per il run :
! daniele, 20/01/2006
! Inserisco flag INOPERT che è settato a 1 per saltare perturbazioni
       write (6,*) 'inopert = ',inopert
       call flush(6)
       IF (INOPERT.EQ.1) GO TO 116
       IM0=1
       IM1=MY
!       IM0=0
!       IM1=0
c se sim 2d inizializzato solo il modo m=m2d
       IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
          IM0=M2D
          IM1=M2D
       ENDIF
! Daniele, 20/02/2009
! per perturbare solo alcuni modi
!       IM0=1
!       IM1=1
       DO 113 IM=IM0,IM1
             M=MM(IM)
             IN0=NANF(M)
             IN1=NANF(M)+NZ(M)-1
             IF (I2D.EQ.1 .OR. I2D.EQ.2) THEN
                IN0=N2D
                IN1=N2D
             ENDIF
! Daniele, 20/02/2009
! per perturbare solo alcuni modi
!             IN0=-10
!             IN1=-9
!             IN0=-8
!             IN1=-8
             DO 113 N=IN0,IN1
                WRITE(6,*) 'perturbo il modo',M,N
                call flush(6)

                J = JANZ(M,N)
ce  aggiunta:
                J1 = JANZ(-M,-N)
                IF(J .EQ. J0) GO TO 113
! daniele,10/12/2008
! tengo solo la pert in vr, come in pixie
! e la metto come in pixie
                VR(0,J) = (0.,0.)
                DO 36 IX=1,LX
                   VR(IX,J) =  0.5*EPS1/X(IX)*sin(pi*X(IX))**3.
! Daniele, 06/05/2009
! aggiunta per avere modi con ampiezze diverse
!                   VR(IX,J) =  VR(IX,J)**((-N-8)**-0.5)
!                   VR(IX,J) =  VR(IX,J)**((-N-8)**-0.25)

!                DO 36 IX=0,LX
!                   VR(IX,J) =  EPS1*(1.0-X(IX)**2)**2
!                   IF(M .EQ. 0) VR(IX,J) = VR(IX,J) * X(IX)
!                   IF(M .GT. 1) VR(IX,J) = VR(IX,J) * X(IX)**(M-1)
                   VR(IX,J1) = CONJG(VR(IX,J))
 36             CONTINUE
! daniele,10/12/2008
! tengo solo la pert in vr, come in pixie
! commentando il seguito
c$$$                IF(M.GT.0)THEN
c$$$                   DO 114 IX=0,LX-1
c$$$                      VT(IX,J)=I*(X(IX+1)*VR(IX+1,J)
c$$$     .                       -X(IX)*VR(IX,J))*ZS2(IX)*2.0/M
c$$$                      VT(IX,J1) = CONJG(VT(IX,J))
c$$$ 114               CONTINUE
c$$$                   VT(-1,J) = VT(0,J)
c$$$                   VT(-1,J1) =CONJG(VT(-1,J))
c$$$                   VT(LX,J) = VT(LX-1,J)*X2(LX-1)/(2.0-X2(LX-1))
c$$$                   VT(LX,J1) = CONJG(VT(LX,J))
c$$$                ELSE
c$$$                   DO 115 IX=0,LX-1
c$$$                      VZ(IX,J)=I*(X(IX+1)*VR(IX+1,J)
c$$$     .                       -X(IX)*VR(IX,J))*ZS2(IX)*2.0*RR/X2(IX)/N
c$$$                      VZ(IX,J1) = CONJG(VZ(IX,J))               
c$$$ 115               CONTINUE                                      
c$$$                   VZ(-1,J) = VZ(0,J)
c$$$                   VZ(-1,J1) =CONJG(VZ(-1,J))
c$$$                   VZ(LX,J) = VZ(LX-1,J)
c$$$                   VZ(LX,J1) = CONJG(VZ(LX,J))
c$$$                ENDIF
 113         CONTINUE
 116         CONTINUE

! daniele, 10/12/2008
! scrivo anche l'istante t=0 che contiene l'equilibrio più la perturbazione assegnata
! oppure se stavo partendo da un'altra simulazione, ho la configurazione da cui parto
C........................................ WRITE ON DISK                 
       WRITE(IBOUT1)  ZEIT,LX,IZ                                        
       DO 1252 J=0,IZ                                                    
          AVT = AIMAG(VT(-1,J))                                         
          AVZ = AIMAG(VZ(-1,J))                                         
          ABT =  REAL(BT(-1,J))                                         
          ABZ =  REAL(BZ(-1,J))                                         
                                                                        
          WRITE(IBOUT1)     AVT, AVZ, ABT, ABZ                          
                                                                        
       DO 1252 IX=0,LX                                                   
          AVR =  REAL(VR(IX,J))                                         
          AVT = AIMAG(VT(IX,J))                                         
          AVZ = AIMAG(VZ(IX,J))                                         
          ABR = AIMAG(BR(IX,J))                                         
          ABT =  REAL(BT(IX,J))                                         
          ABZ =  REAL(BZ(IX,J))                                         
C         ADE =  REAL(DE(IX,J))                                         
C         APR =  REAL(PR(IX,J))                                         
                                                                        
C         IF ( IT .EQ. ITEND ) WRITE(IBDE )   ADE                       
C         IF ( IT .EQ. ITEND ) WRITE(IBPR )   APR                       
C         WRITE(IBDE1)   ADE                                            
C         WRITE(IBPR1)   APR                                            
                                                                        
          WRITE(IBOUT1) AVR, AVT, AVZ, ABR, ABT, ABZ                    
                                                                        
 1252   CONTINUE                                                         

       call flush(IBOUT1)


C.................. MAIN LOOP STARTS HERE...............................
                                                                        
C       WRITE(6,*)'INIZIO MAIN LOOP'                                    

c..... bza definizione:
             bzaini = BZ(LX,0)
             if (ippcd.ne.0) then
!       bzafin = bzaini
                write(6,*) ' ippcd = ',ippcd
                write(6,*) ' bza_ini = ',bzaini
                write(6,*) ' bza_fin = ',bzafin
                call flush(6)                
             endif
                                                                        
       DO 50 IT=1,ITEND                                                 
          ZEIT = ZEIT+DT                                                   

! daniele, 16/07/2008
! per variare la dissipazione durante il ciclo opcd
!          if (ippcd.eq.2) then
!bzn(lx,j) =
!     .          bzaini - ( bzaini-bzafin )*sin(4.*pi*it/itend)
!          endif

        acontr = (bt(50,j0) * conjg(bt(50,j0)))
        afluct = (bt(50,1) * conjg(bt(50,1)))
        if(acontr.gt.1.e5) write(6,*)'acontr fuori contr ! it=',it
        if(acontr.gt.1.e5) go to 5000

!       WRITE(6,*) '   0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0))) = ',
!     .        0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0)))

!          write(6,*) ' it = ',it
!     .         , ', acontr = ',acontr
!     .         ,' afluct = ',afluct
!          call flush(6)      
                                                                        
!       IF (ILINEAR.EQ.1) THEN
! daniele, 28/12/2006
! modifica per calcolo lineare da equilibrio ohmico
       IF (ILINEAR.EQ.1.AND..NOT.IBAND) THEN
         DO 783 IX=0,LX                                                 
!              VR(IX,J0) = (0.0,0.0)                                     
              BR(IX,J0) = (0.0,0.0)                                     
 783     CONTINUE                                                       
         DO 784 IX=-1,LX                                                
!              VT(IX,J0) = (0.0,0.0)                                     
!              VZ(IX,J0) = (0.0,0.0)                                     
              BT(IX,J0) = BT0(IX)                                       
              BZ(IX,J0) = BZ0(IX)                                       
784      CONTINUE                                                       
       END IF                                                           

! daniele, 28/12/2006
! per calcolo lineare da equilibrio ohmico
       IF (ILINEAR.EQ.1.AND.IBAND) THEN
          DO IX=0,LX                                                 
             VR(IX,J0) = vrEQ(ix)
             BR(IX,J0) = brEQ(ix)*I
          ENDDO
          DO IX=-1,LX                                                
             VT(IX,J0) = vtEQ(ix)*I
             VZ(IX,J0) = vzEQ(ix)*I
             BT(IX,J0) = btEQ(ix)
             BZ(IX,J0) = bzEQ(ix)
          ENDDO
       END IF                                                           

! daniele, 11/10/2007
! per tenere costanti le autofunzioni 10, 21, 32 e 43
       J10 = JANZ(1,0)
       J20 = JANZ(2,0)
!       write(6,*) 'br(10,j10) = ',br(10,j10)
!       write(6,*) 'bz(lx,j10) = ',bz(lx,j10)
       IF (ISTATIC10.EQ.1.AND.AF10.NE.0.) THEN
          DO IX=0,LX                                                 
             VR(IX,J10) = vr10(ix)
             BR(IX,J10) = br10(ix)*I
             VR(IX,-J10) = CONJG(VR(IX,J10))
             BR(IX,-J10) = CONJG(BR(IX,J10))
!             VR(IX,J20) = vr20(ix)
!             BR(IX,J20) = br20(ix)*I
!             VR(IX,-J20) = CONJG(VR(IX,J20))
!             BR(IX,-J20) = CONJG(BR(IX,J20))
          ENDDO
          DO IX=-1,LX                                                
             VT(IX,J10) = vt10(ix)*I
             VZ(IX,J10) = vz10(ix)*I
             BT(IX,J10) = bt10(ix)
             BZ(IX,J10) = bz10(ix)
             VT(IX,-J10) = CONJG(VT(IX,J10))
             VZ(IX,-J10) = CONJG(VZ(IX,J10))
             BT(IX,-J10) = CONJG(BT(IX,J10))
             BZ(IX,-J10) = CONJG(BZ(IX,J10))
!             VT(IX,J20) = vt20(ix)*I
!             VZ(IX,J20) = vz20(ix)*I
!             BT(IX,J20) = bt20(ix)
!             BZ(IX,J20) = bz20(ix)
!             VT(IX,-J20) = CONJG(VT(IX,J20))
!             VZ(IX,-J20) = CONJG(VZ(IX,J20))
!             BT(IX,-J20) = CONJG(BT(IX,J20))
!             BZ(IX,-J20) = CONJG(BZ(IX,J20))
          ENDDO
       END IF                                                           

!       IF (ISTATIC10.EQ.1.AND.AF10.EQ.0..AND.IT.EQ.1) THEN
       IF (ISTATIC10.EQ.1.AND.AF10.EQ.0.) THEN
          DO IX=-1,LX                                                
             BZ(IX,J10) = bz10(ix)
             BZ(IX,-J10) = CONJG(BZ(IX,J10))
!             write(6,*) 'bz(ix,j10) = ',ix,bz(ix,j10)
          ENDDO
       END IF                                                           

       J21 = JANZ(2,-1)
       IF (ISTATIC21.EQ.1) THEN
          DO IX=0,LX                                                 
             VR(IX,J21) = vr21(ix)
             BR(IX,J21) = br21(ix)*I
             VR(IX,-J21) = CONJG(VR(IX,J21))
             BR(IX,-J21) = CONJG(BR(IX,J21))
          ENDDO
          DO IX=-1,LX                                                
             VT(IX,J21) = vt21(ix)*I
             VZ(IX,J21) = vz21(ix)*I
             BT(IX,J21) = bt21(ix)
             BZ(IX,J21) = bz21(ix)
             VT(IX,-J21) = CONJG(VT(IX,J21))
             VZ(IX,-J21) = CONJG(VZ(IX,J21))
             BT(IX,-J21) = CONJG(BT(IX,J21))
             BZ(IX,-J21) = CONJG(BZ(IX,J21))
          ENDDO
       END IF                                                           

       J32 = JANZ(3,-2)
       IF (ISTATIC32.EQ.1) THEN
          DO IX=0,LX                                                 
             VR(IX,J32) = vr32(ix)
             BR(IX,J32) = br32(ix)*I
             VR(IX,-J32) = CONJG(VR(IX,J32))
             BR(IX,-J32) = CONJG(BR(IX,J32))
          ENDDO
          DO IX=-1,LX                                                
             VT(IX,J32) = vt32(ix)*I
             VZ(IX,J32) = vz32(ix)*I
             BT(IX,J32) = bt32(ix)
             BZ(IX,J32) = bz32(ix)
             VT(IX,-J32) = CONJG(VT(IX,J32))
             VZ(IX,-J32) = CONJG(VZ(IX,J32))
             BT(IX,-J32) = CONJG(BT(IX,J32))
             BZ(IX,-J32) = CONJG(BZ(IX,J32))
          ENDDO
       END IF                                                           

       J43 = JANZ(4,-3)
       IF (ISTATIC43.EQ.1) THEN
          DO IX=0,LX                                                 
             VR(IX,J43) = vr43(ix)
             BR(IX,J43) = br43(ix)*I
             VR(IX,-J43) = CONJG(VR(IX,J43))
             BR(IX,-J43) = CONJG(BR(IX,J43))
          ENDDO
          DO IX=-1,LX                                                
             VT(IX,J43) = vt43(ix)*I
             VZ(IX,J43) = vz43(ix)*I
             BT(IX,J43) = bt43(ix)
             BZ(IX,J43) = bz43(ix)
             VT(IX,-J43) = CONJG(VT(IX,J43))
             VZ(IX,-J43) = CONJG(VZ(IX,J43))
             BT(IX,-J43) = CONJG(BT(IX,J43))
             BZ(IX,-J43) = CONJG(BZ(IX,J43))
          ENDDO
       END IF                                                           
                                                                        
       J = JANZ(MPLO(2),NPLO(2))                                        
       VRALT = VR(IG,J)                                                 
       VZALT = VZ(IG,J)                                                 
C      VZALT = 1.0                                                      
C      VRALT = 1.0                                                      
                                                                        
C............... B UND P - HALBSCHRITT .................................
                                                                        
       DO 12 IM=0,MY                                                    
       M = MM(IM)                                                       
       DO 12 N=NA1(M),NE(M)                                             
       DO 55 IX=0,LX                                                    
       SUM0(IX) = (0.0,0.0)                                             
       SUM1(IX) = (0.0,0.0)                                             
       SUM2(IX) = (0.0,0.0)                                             
       SUM3(IX) = (0.0,0.0)                                             
       SUM4(IX) = (0.0,0.0)                                             
       SUM5(IX) = (0.0,0.0)                                             
       SUM6(IX) = (0.0,0.0)                                             
       SUM7(IX) = (0.0,0.0)                                             
       SUM8(IX) = (0.0,0.0)                                             
       SUM9(IX) = (0.0,0.0)                                             
       SUM10(IX) = (0.0,0.0)                                            
   55 CONTINUE                                                          
       DO 56 MS=-MM(MY)+M,MM(MY)                                        
       IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) GO TO 56                  
       NS0 = NA(MS)
       NS1 = NE(MS)
       IF (I2D.EQ.1) THEN
          NS0 = MS*HELI2D
          NS1 = MS*HELI2D  
       ENDIF
       DO 57 NS=NS0,NS1
       IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) GO TO 57
       J = JANZ(MS,NS)                                                  
       J1 = JANZ(M-MS,N-NS)                                             
!       write (6,*) 'm,n,ms,ns,m-ms,n-ns',m,n,ms,ns,m-ms,n-ns
!       call flush(6)
       DO 58 IX=0,LX-1                                                  
       SUM3(IX) = SUM3(IX)+VT(IX,J)*BZ(IX,J1)-VZ(IX,J)*BT(IX,J1)        
   58 CONTINUE                                                          
       DO 26 IX=0,LX                                                    
       SUM0(IX) = SUM0(IX)+BR(IX,J)*(VT(IX-1,J1)+VT(IX,J1))             
       SUM1(IX) = SUM1(IX)+VR(IX,J)*(BT(IX-1,J1)+BT(IX,J1))             
       SUM2(IX) = SUM2(IX)+BR(IX,J)*(VZ(IX-1,J1)+VZ(IX,J1))             
       SUM8(IX) = SUM8(IX)+VR(IX,J)*(BZ(IX-1,J1)+BZ(IX,J1))             
   26 CONTINUE                                                          
   57 CONTINUE                                                          
   56 CONTINUE                                                          
C ACHTUNG MANCHE SUMMEN (ABL.) AUCH AM RAND NOTWENDIG                   
       DO 59 IX=0,LX-1                                                  
       DSUM2(IX)=(SUM1(IX+1)-SUM1(IX)-SUM0(IX+1)+SUM0(IX))*ZS22(IX)*0.5 
       DSUM3(IX) =((SUM2(IX+1)-SUM8(IX+1))*X(IX+1)-(SUM2(IX)-SUM8(IX))  
     1 *X(IX))*ZS22(IX)/X2(IX)*0.5                                      
   59 CONTINUE                                                          
       J = JANZ(MM(IM),N)                                               
       J1 = JANZ(-MM(IM),-N)                                            
       CFAK1 = N/RR*I                                                   
       DO 14 IX=0,LX-1                                                  
       CFAK = M/X2(IX)*I                                                
       BTN(IX,J) = BT(IX,J)+ADT*(+CFAK1*SUM3(IX)-DSUM2(IX))             
       BZN(IX,J) = BZ(IX,J)+ADT*(+DSUM3(IX)-CFAK*SUM3(IX))              
   14 CONTINUE                                                          
!       BZN(LX,J) = BZN(LX-1,J)
! Daniele, 28/11/2008
! metto le B.C.s nuove anche per i modi diversi da zero
!         BZN(LX,J) = BZN(LX-1,J)*
!     .        (2.*ETAA(LX)+VRLX/LX)/
!     .        (2.*ETAA(LX)-VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
         BZN(LX,J) = BZN(LX-1,J)*
     .        (2.*ETAA(LX)+VRLX/LX)/
     .        (2.*ETAA(LX)-VRLX/LX)
     .        + I*ETAA(LX)*N/RR*BR(LX,J)*
     .        (2./LX)/(2.*ETAA(LX)-VRLX/LX)

!       BTN(LX,J) = BTN(LX-1,J)*X2(LX-1)/(2.-X2(LX-1))                   
! Daniele, 28/11/2008
! metto le B.C.s nuove anche per i modi diversi da zero
!       BTN(LX,J) = BTN(LX-1,J)*
!     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
!     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
       BTN(LX,J) = BTN(LX-1,J)*  
     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
     .        + I*ETAA(LX)*M*BR(LX,J)*
     .        (2./LX)/(2.*X2(LX)*ETAA(LX)-VRLX/LX)

       BTN(-1,J) = -BTN(0,J)                                            
       BZN(-1,J) = -BZN(0,J)                                            
       IF(J.NE.0) GO TO 66                                              
!         BTN(LX,J) = BT0(LX)
! Daniele, 10/05/2007
!         BTN(LX,J) = 2.*BTWALL-BTN(LX-1,J)
! Daniele, 25/11/2008
! invece di mettere condizione per Ip costante, metto E0 costante
!         BTN(LX,J) = BTN(LX-1,J)*X2(LX-1)/X2(LX)
!     .        + E0/(ETAA(LX)*LX*X2(LX))
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
       BTN(LX,J) = BTN(LX-1,J)*
     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
     .        + E0*(2./LX)/(2.*X2(LX)*ETAA(LX)-VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
! nel caso J=0 non cambia nulla!

C     ............B.C. PER BZ ALLA PARETE:                              
!         BZN(LX,J) = BZN(LX-1,J)
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
         BZN(LX,J) = BZN(LX-1,J)*
     .        (2.*ETAA(LX)+VRLX/LX)/
     .        (2.*ETAA(LX)-VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
! nel caso J=0 non cambia nulla!

! daniele, 11/05/2007
! condizioni al bordo per vr finito
!         IF (IPR.EQ.-2) THEN
!            BZN(LX,J) = (2.*ETA0+ZDR*VR0(LX))/(2.*ETA0-ZDR*VR0(LX))*
!     .             BZN(LX-1,J)
!         END IF
C        BZN(LX,J) = BZ0(LX)                                            
         if (ippcd.eq.1) bzn(lx,j) =
     .          bzaini - it * ( bzaini-bzafin )/itend
         if (ippcd.eq.2) bzn(lx,j) =
     .          bzaini - ( bzaini-bzafin )*sin(4.*pi*it/itend)
         if (ippcd.eq.3) bzn(lx,j) = bzaini
         if (ippcd.eq.4) bzn(lx,j) =
     .          bzaini - ( bzaini-bzafin )*
     .          (1.-cos(4.*pi*it/itend))/2.

!       WRITE(6,*) 'bc1   0.5*(REAL(BTN(LX-1,J0))+REAL(BTN(LX,J0))) = ',
!     .        0.5*(REAL(BTN(LX-1,J0))+REAL(BTN(LX,J0)))

C    ........                                                           
66       IF(M .NE. 0) GO TO 68                                          
         BZN(-1,J) = BZN(0,J)                                           
68     J1 = JANZ(-M,-N)                                                 
                                                                        
C  BR AUS BT +BZ ........................                               
       BRN(0,J) =(0.0,0.0)                                              
       DO 149 IX=1,LX-1                                                 
       BRN(IX,J)=(I/ZS22(IX-1)*(-BTN(IX-1,J)*M-BZN(IX-1,J)*N*X2(IX-1)   
     1 /RR)+BRN(IX-1,J)*X(IX-1))/X(IX)                                  
  149 CONTINUE                                                          
      IF(M.NE.1) GO TO 67                                               
C      BTN(0,J) = -BRN(0,J)                                             
       BRN(0,J) = (4.0*BRN(1,J)-BRN(2,J))*0.333333333                   
       BTN(-1,J) = -BTN(0,J)+2.0*I*BRN(0,J)                             
   67  IF(J.EQ.0) GO TO 12                                              
       DO 98 IX=0,LX                                                    
       BRN(IX,J1) = CONJG(BRN(IX,J))                                    
   98 CONTINUE                                                          
       DO 401 IX=-1,LX                                                  
       BTN(IX,J1) = CONJG(BTN(IX,J))                                    
       BZN(IX,J1) = CONJG(BZN(IX,J))                                    
  401 CONTINUE                                                          
   12 CONTINUE                                                          
                                                                        
C      WRITE(6,*)' FATTO B*',IT,'   IT'                                 
C .................... J-BESTIMMUNG ................................... 
                                                                        
       DO 410 IM=0,MY                                                   
       M = MM(IM)                                                       
       DO 410 N=NA1(M),NE(M)                                            
       FAK1 = N/RR                                                      
       J = JANZ(M,N)                                                    
       J1 = JANZ(-M,-N)                                                 
       DO 408 IX=0,LX-1                                                 
       JR(IX,J) = I*(M/X2(IX)*BZN(IX,J)-FAK1*BTN(IX,J))                 
  408 CONTINUE                                                          
       DO 409 IX=1,LX-1                                                 
       JT(IX,J) = I*FAK1*BRN(IX,J)-2.0*ZS2(IX)*(BZN(IX,J)-BZN(IX-1,J))  
       JZ(IX,J) = 2.0*ZS2(IX)/X(IX)*(X2(IX)*BTN(IX,J)-X2(IX-1)          
     1 *BTN(IX-1,J))-I*M/X(IX)*BRN(IX,J)                                
  409 CONTINUE                                                          
       JT(LX,J) = (0.0,0.0)                                             
       JZ(LX,J) = (0.0,0.0)                                             
       JT(0,J) =(0.0,0.0)                                               
       JZ(0,J) =(0.0,0.0)                                               
       IF(M.NE.0) GO TO 434                                             
       JT(LX,J) = -2.*ZS2(LX)*(BZN(LX,J) - BZN(LX-1,J))                 
       JZ(LX,J) =  2.*ZS2(LX)*((2. - X2(LX-1))*BTN(LX,J) -              
     &                                          X2(LX-1)*BTN(LX-1,J))   
        JZ(0,J) = 2.*BTN(0,J)/X2(0)                                     
  434 IF(M.NE.1) GO TO 420                                              
       JT(0,J) = JT(1,J)                                                
  420 IF(J.EQ.0) GO TO 410                                              
       DO 437 IX=0,LX                                                   
       JR(IX,J1) = CONJG(JR(IX,J))                                      
  437 CONTINUE                                                          
       DO 421 IX=0,LX                                                   
       JT(IX,J1) = CONJG(JT(IX,J))                                      
       JZ(IX,J1) = CONJG(JZ(IX,J))                                      
  421 CONTINUE                                                          
  410 CONTINUE                                                          
                                                                        
C ............................DENSITY ..................................
                                                                        
C     DO 700 IM=0,MY                                                    
C      M = MM(IM)                                                       
C      DO 700 N=NA1(M),NE(M)                                            
C                                                                       
C       DO 701 IX = 0,LX                                                
C       SUM0(IX) = (0.0,0.0)                                            
C       SUM1(IX) = (0.0,0.0)                                            
C701     CONTINUE                                                       
C                                                                       
C       DO 702 MS = -MM(MY)+M,MM(MY)                                    
C      IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) GO TO 702                 
C        DO 703 NS = NA(MS),NE(MS)                                      
C      IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) GO TO 703               
C        J = JANZ(MS,NS)                                                
C        J1= JANZ(M-MS,N-NS)                                            
C          DO 704 IX=1,LX-1                                             
C                                                                       
C          VTIXJ= (VT(IX-1,J)+VT(IX,J))/2.                              
C          VZIXJ= (VZ(IX-1,J)+VZ(IX,J))/2.                              
C                                                                       
C     DENSITY :                                                         
C          S0 = - ZDR * 0.5 / X(IX)                                     
C          S1 = X(IX+1)*DE(IX+1,J1)*VR (IX+1,J)                         
C          S2 = X(IX-1)*DE(IX-1,J1)*VR (IX-1,J)                         
C          S3 = DE(IX,J1)                                               
C          S4 = I * M * VTIXJ  / X(IX)                                  
C          S5 = I * N * VZIXJ   / RR                                    
C          SUM0(IX) = S0 *(S1-S2) - S3 *(S4+S5) + SUM0(IX)              
C                                                                       
C    PRESSURE :                                                         
C          S0 = - ZDR * 0.5 / X(IX)                                     
C          S1 = X(IX+1)*PR(IX+1,J1)*VR (IX+1,J)                         
C          S2 = X(IX-1)*PR(IX-1,J1)*VR (IX-1,J)                         
C          S3 = PR(IX,J1)                                               
C          S4 = I * M * VTIXJ  / X(IX)                                  
C          S5 = I * N * VZIXJ   / RR                                    
C          S6 = -2./3. * PR(IX,J1) * (                                  
C    .    1./X(IX) * (X(IX+1)*VR(IX+1,J)-X(IX-1)*VR(IX-1,J)) *0.5 *ZDR  
C    .       + I * MS /X(IX) * VTIXJ + I * NS /RR * VZIXJ    )          
C          SUM1(IX) = S0 *(S1-S2) - S3 *(S4+S5) + S6 + SUM1(IX)         
C                                                                       
C704       CONTINUE                                                     
C703       CONTINUE                                                     
C702       CONTINUE                                                     
C                                                                       
C         J = JANZ(MM(IM),N)                                            
C         J1= JANZ(-MM(IM),-N)                                          
C         DO 705 IX=1,LX-1                                              
C                                                                       
C          AMX = 1.0                                                    
C          DENN(IX,J)=AMX*(DE(IX+1,J)+DE(IX-1,J))/2.+ (1.-AMX)*DE(IX,J) 
C    .                + ADT * SUM0(IX)                                  
C                                                                       
C          PRNN(IX,J)=AMX*(PR(IX+1,J)+PR(IX-1,J))/2.+ (1.-AMX)*PR(IX,J) 
C    .                + ADT * SUM1(IX)                                  
C                                                                       
C          DENN(IX,J1)= CONJG(DENN(IX,J))                               
C          PRNN(IX,J1)= CONJG(PRNN(IX,J))                               
C705       CONTINUE                                                     
C          IF ( MM(IM) .EQ. 0 ) THEN                                    
C          DENN(0,J)  = DENN(1,J)                                       
C          DENN(0,J1) = CONJG(DENN(0,J))                                
C          PRNN(0,J)  = PRNN(1,J)                                       
C         PRNN(0,J1) = CONJG(PRNN(0,J))                                 
C         ELSE                                                          
C         DENN(0,J)  = 0.0                                              
C         DENN(0,J1) = CONJG(DENN(0,J))                                 
C         PRNN(0,J)  = 0.0                                              
C         PRNN(0,J1) = CONJG(PRNN(0,J))                                 
C         ENDIF                                                         
C         IF( J .EQ. 0 )THEN                                            
C         DENN(LX,J) = DENN(LX-1,J)                                     
C         ELSE                                                          
C         DENN(LX,J) = 0.0                                              
C         ENDIF                                                         
C         DENN(LX,J1)=CONJG(DENN(LX,J))                                 
C         PRNN(LX,J) = PRNN(LX-1,J)                                     
C         PRNN(LX,J) = PR(LX,J)                                         
C         PRNN(LX,J1)=CONJG(PRNN(LX,J))                                 
C700       CONTINUE                                                     
C                                                                       
C ............  GANZSCHRITTSEMI-IMPLIZ. VERFAHREN ......................
                                                                        
       DO 15 IM=0,MY                                                    
       M = MM(IM)                                                       
       DO 15 N=NA1(M),NE(M)                                             
       ampv = 1.e0
c       if(m.eq.1 .and. n.eq.-10) ampv= 30.e0/2.e0       
c       write(6,*) m,n
c       if(m.eq.1 .and. n.eq.-10) write(6,*) 'ampv= ',ampv  
       FAK1 = N/RR                                                      
       DO 16 IX=0,LX                                                    
       SUM0(IX) = (0.0,0.0)                                             
       SUM1(IX) = (0.0,0.0)                                             
       SUM2(IX) = (0.0,0.0)                                             
       SUM3(IX) = (0.0,0.0)                                             
       SUM4(IX) = (0.0,0.0)                                             
       SUM5(IX) = (0.0,0.0)                                             
       SUM6(IX) = (0.0,0.0)                                             
       SUM7(IX) = (0.0,0.0)                                             
       SUM8(IX) = (0.0,0.0)                                             
       SUM9(IX) = (0.0,0.0)                                             
       SUM10(IX) = (0.0,0.0)                                            
       SUM11(IX) = (0.0,0.0)                                            
       SUM12(IX) = (0.0,0.0)                                            
       SUM13(IX) = (0.0,0.0)                                            
       SUM14(IX) = (0.0,0.0)                                            
       SUM15(IX) = (0.0,0.0)                                            
       SUM16(IX) = (0.0,0.0)                                            
! daniele, 24/11/2006
! aggiungo per escludere il termine convettivo
       DSUM1(IX) = (0.0,0.0)
       DSUM2(IX) = (0.0,0.0)
       DSUM3(IX) = (0.0,0.0)
   16 CONTINUE                                                          
       DO 17 MS=-MM(MY)+M,MM(MY)                                        
       IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) GO TO 17                  
       NS0 = NA(MS)
       NS1 = NE(MS)
       IF (I2D.EQ.1) THEN
          NS0 = MS*HELI2D
          NS1 = MS*HELI2D  
       ENDIF
       DO 18 NS=NS0,NS1
       IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) GO TO 18
       J = JANZ(MS,NS)                                                  
       J1 = JANZ(M-MS,N-NS)                                             
       F1 = NS/RR                                                       
       DO 19 IX=1,LX-1                                                  
       SUM1(IX) = SUM1(IX)+(VR(IX+1,J)-VR(IX-1,J))*VR(IX,J1)            
       SUM2(IX)=SUM2(IX)+(VT(IX,J)+VT(IX-1,J))*(VT(IX,J1)+VT(IX-1,J1))  
       SUM3(IX) = SUM3(IX)+MS*VR(IX,J)*(VT(IX,J1)+VT(IX-1,J1))          
       SUM4(IX) = SUM4(IX)+NS*VR(IX,J)*(VZ(IX,J1)+VZ(IX-1,J1))          
       SUM12(IX) = SUM12(IX)+JT(IX,J)*(BZN(IX,J1)+BZN(IX-1,J1))         
       SUM0(IX) = SUM0(IX)+JZ(IX,J)*(BTN(IX,J1)+BTN(IX-1,J1))           
   19 CONTINUE                                                          
       DO 27 IX=1,LX-2                                                  
       SUM5(IX)=SUM5(IX)+(VT(IX+1,J)-VT(IX-1,J))*(VR(IX,J1)+VR(IX+1,J1))
       SUM6(IX)=SUM6(IX)+VT(IX,J)*(VR(IX,J1)+VR(IX+1,J1))               
       SUM7(IX)=SUM7(IX)+MS*VT(IX,J)*VT(IX,J1)                          
       SUM8(IX)=SUM8(IX)+NS*VT(IX,J)*VZ(IX,J1)                          
       SUM9(IX)=SUM9(IX)+(VZ(IX+1,J)-VZ(IX-1,J))*(VR(IX,J1)+VR(IX+1,J1))
       SUM10(IX)=SUM10(IX)+MS*VZ(IX,J)*VT(IX,J1)                        
       SUM11(IX)=SUM11(IX)+NS*VZ(IX,J)*VZ(IX,J1)                        
   27 CONTINUE                                                          
       DO 423 IX=0,LX-1                                                 
       SUM13(IX) = SUM13(IX)+(JZ(IX+1,J)+JZ(IX,J))*(BRN(IX+1,J1)+       
     1 BRN(IX,J1))                                                      
       SUM14(IX) = SUM14(IX)+JR(IX,J)*BZN(IX,J1)                        
       SUM15(IX) = SUM15(IX)+JR(IX,J)*BTN(IX,J1)                        
       SUM16(IX) = SUM16(IX)+(JT(IX+1,J)+JT(IX,J))*(BRN(IX+1,J1)+       
     1 BRN(IX,J1))                                                      
  423 CONTINUE                                                          
   18 CONTINUE                                                          
   17 CONTINUE                                                          
! daniele, 24/11/2006
! aggiungo per escludere il termine convettivo
      IF (INOCONV.EQ.1) GO TO 399
       DO 429 IX=1,LX-1                                                 
      DSUM1(IX) = ZS2(IX)*SUM1(IX)-0.25/X(IX)*SUM2(IX)+0.5*I*(SUM3(IX)  
     3 /X(IX)+SUM4(IX)/RR)                                              
  429 CONTINUE                                                          
       DO 424 IX=1,LX-2                                                 
       DSUM2(IX)=(SUM5(IX)*ZS22(IX)*0.5+SUM6(IX)/X2(IX))*0.5+I*(SUM7(IX)
     1 /X2(IX)+SUM8(IX)/RR)                                             
       DSUM3(IX) = SUM9(IX)*ZS22(IX)*0.25+I*(SUM10(IX)/                 
     2 X2(IX)+SUM11(IX)/RR)                                             
  424 CONTINUE                                                          
 399  DSUM2(LX-1) = DSUM2(LX-2)                                        
       DSUM3(LX-1) = DSUM3(LX-2)                                        
       DSUM2(0) =(0.0,0.0)                                              
       DSUM3(0) =(0.0,0.0)                                              
       IF(M.NE.1) GO TO 411                                             
       DSUM2(0) = DSUM2(1)                                              
  411  IF(M.NE.0) GO TO 412                                             
       DSUM3(0) = DSUM3(1)                                              
  412  J = JANZ(M,N)                                                    
       J1 = JANZ(-M,-N)                                                 
                                                                        
       DO 20 IX=1,LX-1                                                  
C        RHOX = REAL (DENN(IX,0))                                       
C        IF( IACC .NE. 1 ) RHOX = 1.                                    
         KMN1 = (M/X(IX))**2 + (N/RR)**2                                
         R1(IX) = VR(IX,J)                                              
     1   - A22*( ZDR**2*(VR(IX+1,J)-2.*VR(IX,J)+VR(IX-1,J))             
     2        + ZDR*0.5/X(IX)*(VR(IX+1,J)-VR(IX-1,J))                   
     3        - (KMN1+1./ X(IX)**2)*VR(IX,J)                            
     4        - I*M/X(IX)**2*(VT(IX,J)+VT(IX-1,J)) )                    
     5   - DT*( DSUM1(IX) +  0.5*(-SUM12(IX) + SUM0(IX)) )              
C    5   - DT*( DSUM1(IX) +  0.5*(-SUM12(IX) + SUM0(IX)                 
C    6          + BPRES * (PRNN(IX+1,J)-PRNN(IX-1,J)) * ZDR             
C    7           )/RHOX )                                               
   20 CONTINUE                                                          
                                                                        
       DO 430 IX=0,LX-1                                                 
C        RHOX2 = ( REAL(DENN(IX,0))+REAL(DENN(IX+1,0)) )/2.             
C        IF ( IACC .NE. 1 ) RHOX2 = 1.                                  
         KMN2 = (M/X2(IX))**2 + (N/RR)**2                               
         R2(IX) = VT(IX,J)                                              
     1   - A22 * ( ZDR**2/X2(IX)*                                       
     2     (X(IX+1)*(VT(IX+1,J)-VT(IX,J))-X(IX)*(VT(IX,J)-VT(IX-1,J)))  
     3     - (KMN2+1./X2(IX)**2)*VT(IX,J)                               
     4     + I*M/X2(IX)**2*(VR(IX+1,J)+VR(IX,J))  )                     
     5  - DT * ( DSUM2(IX)  -SUM13(IX)*0.25+SUM14(IX)  )                
C    5  - DT * ( DSUM2(IX) +( -SUM13(IX)*0.25+SUM14(IX)                 
C    6    + BPRES * I * M /X2(IX) *(PRNN(IX-1,J)+PRNN(IX+1,J))/2.       
C    7           )/RHOX2 )                                              
         R3(IX) = VZ(IX,J)                                              
     1   - A22 * ( ZDR**2/X2(IX)*                                       
     2      (X(IX+1)*(VZ(IX+1,J)-VZ(IX,J))-X(IX)*(VZ(IX,J)-VZ(IX-1,J))) 
     3      - KMN2*VZ(IX,J)  )                                          
     4   - DT * ( DSUM3(IX)  -SUM15(IX)+SUM16(IX)*0.25 )                
C    4   - DT * ( DSUM3(IX) +( -SUM15(IX)+SUM16(IX)*0.25                
C    5      + BPRES * I * N / RR *(PRNN(IX-1,J)+PRNN(IX+1,J))/2.        
C    6           )/RHOX2 )                                              
  430 CONTINUE                                                          
                                                                        
         DO 801 L=1,2                                                   
              FF(L,0) = (0.0,0.0)                                       
         DO 801 K=1,2                                                   
              EE(L,K,0) = (0.0,0.0)                                     
 801     CONTINUE                                                       
                                                                        
         A12 = ampv * MU2(0)*DT + A22                                          
C        RHO02 = (REAL(DENN(0,0))+REAL(DENN(1,0)))/2.                   
C        IF( IACC .NE. 1 ) RHO02 = 1.                                   
C        A12 = MU2(0)*DT/RHO02 + A22                                    
         KMN2 = (M/X2(0))**2 + (N/RR)**2                                
         ALFS = 1. + A12*ZDR**2*X(1)/X2(0) + A12*(KMN2+1./X2(0)**2)     
         BETS = 1. + A12*ZDR**2*X(1)/X2(0) + A12*KMN2                   
         EE(2,1,0) = I*M*A12/ALFS/X2(0)**2                              
         EE(2,2,0) = A12*ZDR**2*X(1)/ALFS/X2(0)                         
         IF (M.EQ.1) THEN                                               
              EE(1,1,0) = (1.0,0.0)                                     
              EE(2,1,0) = 2.*EE(2,1,0)                                  
         END IF                                                         
         FF(2,0) = R2(0)/ALFS                                           
         AAA(0) = A12*ZDR**2/BETS*X(1)/X2(0)                            
         BBB(0) = R3(0)/BETS                                            
                                                                        
         DO 802 IX=1,LX-1                                               
              KMN1 = (M/X(IX))**2 + (N/RR)**2                           
              KMN2 = (M/X2(IX))**2 + (N/RR)**2                          
C             RHOX1= REAL(DENN(IX,0))                                   
C             RHOX2= (REAL(DENN(IX,0))+REAL(DENN(IX+1,0)))/2.           
C             IF ( IACC .NE. 1 ) RHOX1= 1.                              
C             IF ( IACC .NE. 1 ) RHOX2= 1.                              
              A11 = ampv * MU1(IX)*DT + A22                                    
              A12 = ampv * MU2(IX)*DT + A22                                    
C             A11 = MU1(IX)*DT/RHOX1 + A22                              
C             A12 = MU2(IX)*DT/RHOX2 + A22                              
                                                                        
           AA(1,1,IX) = -A11*(ZDR**2 + ZDR*0.5/X(IX))                   
           AA(1,2,IX) = (0.0,0.0)                                       
           AA(2,1,IX) = -A12*I*M/X2(IX)**2                              
           AA(2,2,IX) = -A12*ZDR**2*X(IX+1)/X2(IX)                      
                                                                        
           BB(1,1,IX) = 1. + A11*(2.*ZDR**2 + (KMN1 + 1./X(IX)**2))     
           BB(1,2,IX) =  A11*I*M/X(IX)**2                               
           BB(2,1,IX) = -A12*I*M/X2(IX)**2                              
           BB(2,2,IX) = 1. + A12*(2.*ZDR**2 + (KMN2 + 1./X2(IX)**2))    
                                                                        
           CC(1,1,IX) = -A11*(ZDR**2 - ZDR*0.5/X(IX))                   
           CC(1,2,IX) =  A11*I*M/X(IX)**2                               
           CC(2,1,IX) = (0.0,0.0)                                       
           CC(2,2,IX) = -A12*ZDR**2*X(IX)/X2(IX)                        
                                                                        
           ALJ(IX)  = -A12*ZDR**2*X(IX+1)/X2(IX)                        
           BETJ(IX) = 1. + A12*2.*ZDR**2 + A12*KMN2                     
           GAMJ(IX) = -A12*ZDR**2*X(IX)/X2(IX)                          
                                                                        
           DD(1,IX) = R1(IX)                                            
           DD(2,IX) = R2(IX)                                            
                                                                        
              DO 803 L=1,2                                              
                   AV1(L) = FF(L,IX-1)                                  
                   AV2(L) = DD(L,IX)                                    
              DO 803 K=1,2                                              
                   AM1(L,K) = AA(L,K,IX)                                
                   AM2(L,K) = BB(L,K,IX)                                
                   AM3(L,K) = CC(L,K,IX)                                
                   AM4(L,K) = EE(L,K,IX-1)                              
 803          CONTINUE                                                  
              CALL MATCAL(AM1,AM2,AM3,AM4,AM5,AV1,AV2,AV3)              
              DO 804 L=1,2                                              
                   FF(L,IX) = AV3(L)                                    
              DO 804 K=1,2                                              
                   EE(L,K,IX) = AM5(L,K)                                
 804          CONTINUE                                                  
              AAA(IX) = -ALJ(IX)/(BETJ(IX) + GAMJ(IX)*AAA(IX-1))        
              BBB(IX) = (R3(IX) - GAMJ(IX)*BBB(IX-1))/                  
     &                           (BETJ(IX) + GAMJ(IX)*AAA(IX-1))        
 802     CONTINUE                                                       
                                                                        
         UU(1,LX) = (0.0,0.0)                                           

!         UU(2,LX) = (0.0,0.0)                                           
!         VZN(LX,J) = (0.0,0.0)                                          
! daniele, 02/12/2008
! modifiche alle B.C. di vtheta e vz, per averle 0. esattamente a r=1.

! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
         IF (J.EQ.0) UU(1,LX) = VRLX

! daniele, 02/12/2008
! modifiche alle B.C. di vtheta e vz, per averle 0. esattamente a r=1.
         UU(2,LX-1) = (EE(2,1,LX-1)*UU(1,LX) + FF(2,LX-1))
     .          /(1. + EE(2,2,LX-1))
         UU(2,LX) = - UU(2,LX-1)
         UU(1,LX-1) = EE(1,1,LX-1)*UU(1,LX) + EE(1,2,LX-1)*UU(2,LX)
     .          + FF(1,LX-1)

         VZN(LX-1,J) = BBB(LX-1)/(1. + AAA(LX-1))
         VZN(LX,J) = - VZN(LX-1,J)

!         DO 807 IX=LX-1,0,-1
! daniele, 02/12/2008
! modifiche alle B.C. di vtheta e vz, per averle 0. esattamente a r=1.
         DO 807 IX=LX-2,0,-1
              UU(1,IX) = EE(1,1,IX)*UU(1,IX+1) + EE(1,2,IX)*UU(2,IX+1)  
     &                                         + FF(1,IX)               
              UU(2,IX) = EE(2,1,IX)*UU(1,IX+1) + EE(2,2,IX)*UU(2,IX+1)  
     &                                         + FF(2,IX)               
              VZN(IX,J) = AAA(IX)*VZN(IX+1,J) + BBB(IX)                 
 807     CONTINUE                                                       
         DO 808 IX=0,LX                                                 
              VRN(IX,J) = UU(1,IX)                                      
              VTN(IX,J) = UU(2,IX)                                      
 808     CONTINUE                                                       
                                                                        
         VTN(-1,J) = - VTN(0,J)                                         
         VZN(-1,J) = - VZN(0,J)                                         
         IF (M.NE.0) GO TO 809                                          
         VZN(-1,J) = VZN(0,J)                                           
 809     IF (M.NE.1) GO TO 810                                          
         VTN(-1,J) = VTN(0,J)                                           
                                                                        
 810   DO 351 IX=0,LX                                                   
       VRN(IX,J1) = CONJG(VRN(IX,J))                                    
  351 CONTINUE                                                          
       DO 352 IX=-1,LX                                                  
       VTN(IX,J1) = CONJG(VTN(IX,J))                                    
       VZN(IX,J1) = CONJG(VZN(IX,J))                                    
  352 CONTINUE                                                          
 15   CONTINUE                                                          
                                                                        
C ............................  B  .....................................
       DO 28 IM=0,MY                                                    
       M = MM(IM)                                                       
       DO 28 N=NA1(M),NE(M)                                             
       DO 29 IX=0,LX                                                    
       SUM0(IX) = (0.0,0.0)                                             
       SUM1(IX) = (0.0,0.0)                                             
       SUM2(IX) = (0.0,0.0)                                             
       SUM3(IX) = (0.0,0.0)                                             
       SUM8(IX) = (0.0,0.0)                                             
   29 CONTINUE                                                          
       DO 30 MS=-MM(MY)+M,MM(MY)                                        
       IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) GO TO 30                  
       NS0 = NA(MS)
       NS1 = NE(MS)
       IF (I2D.EQ.1) THEN
          NS0 = MS*HELI2D
          NS1 = MS*HELI2D  
       ENDIF
       DO 31 NS=NS0,NS1
       IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) GO TO 31
       J = JANZ(MS,NS)                                                  
       J1 = JANZ(M-MS,N-NS)                                             
       DO 32 IX=0,LX-1                                                  
       SUM3(IX) = SUM3(IX)+(VT(IX,J)+VTN(IX,J))*BZN(IX,J1)              
     1 -(VZ(IX,J)+VZN(IX,J))*BTN(IX,J1)                                 
   32 CONTINUE                                                          
       DO 33 IX=0,LX                                                    
       SUM0(IX) = SUM0(IX)+BRN(IX,J)*                                   
     1 (VT(IX-1,J1)+VTN(IX-1,J1)+VT(IX,J1)+VTN(IX,J1))                  
       SUM1(IX)=SUM1(IX)+(VR(IX,J)+VRN(IX,J))*(BTN(IX-1,J1)+BTN(IX,J1)) 
       SUM2(IX) = SUM2(IX)+BRN(IX,J)*                                   
     1 (VZ(IX-1,J1)+VZN(IX-1,J1)+VZ(IX,J1)+VZN(IX,J1))                  
       SUM8(IX) = SUM8(IX)+(VR(IX,J)+VRN(IX,J))*                        
     1 (BZN(IX-1,J1)+BZN(IX,J1))                                        
   33 CONTINUE                                                          
   31 CONTINUE                                                          
   30 CONTINUE                                                          
C ACHTUNG MANCHE SUMMEN (ABL.) AUCH AM RAND NOTWENDIG                   
       DO 34 IX=0,LX-1                                                  
       DSUM2(IX)=(SUM1(IX+1)-SUM1(IX)-SUM0(IX+1)+SUM0(IX))*ZS22(IX)*0.5 
       DSUM3(IX) =((SUM2(IX+1)-SUM8(IX+1))*X(IX+1)-(SUM2(IX)-SUM8(IX))  
     1 *X(IX))*ZS22(IX)/X2(IX)*0.5                                      
   34 CONTINUE                                                          
       J = JANZ(MM(IM),N)                                               
       J1 = JANZ(-MM(IM),-N)                                            
       CFAK1 = N/RR*I                                                   
       DO 35 IX=0,LX-1                                                  
       CFAK = M/X2(IX)*I                                                
       BT(IX,J) = BT(IX,J)+DT2*(+CFAK1*SUM3(IX)-DSUM2(IX))              
       BZ(IX,J) = BZ(IX,J)+DT2*(+DSUM3(IX)-CFAK*SUM3(IX))               
   35 CONTINUE                                                          
!       BZ(LX,J) = BZ(LX-1,J)
! Daniele, 28/11/2008
! metto le B.C.s nuove anche per i modi diversi da zero
!         BZ(LX,J) = BZ(LX-1,J)*
!     .        (2.*ETAA(LX)+VRLX/LX)/
!     .        (2.*ETAA(LX)-VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
         BZ(LX,J) = BZ(LX-1,J)*
     .        (2.*ETAA(LX)+VRLX/LX)/
     .        (2.*ETAA(LX)-VRLX/LX)
     .        + I*ETAA(LX)*N/RR*BRN(LX,J)*
     .        (2./LX)/(2.*ETAA(LX)-VRLX/LX)

!       BT(LX,J) = BT(LX-1,J)*X2(LX-1)/(2.-X2(LX-1))                   
! Daniele, 28/11/2008
! metto le B.C.s nuove anche per i modi diversi da zero
!       BT(LX,J) = BT(LX-1,J)*
!     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
!     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
       BT(LX,J) = BT(LX-1,J)*
     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
     .        + I*ETAA(LX)*M*BRN(LX,J)*
     .        (2./LX)/(2.*X2(LX)*ETAA(LX)-VRLX/LX)

       BT(-1,J) = -BT(0,J)                                              
       BZ(-1,J) = -BZ(0,J)                                              
       IF(J.NE.0) GO TO 72                                              
!         BT(LX,J) = BT0(LX)
! Daniele, 10/05/2007
!         BT(LX,J) = 2.*BTWALL-BT(LX-1,J)
! Daniele, 25/11/2008
! invece di mettere condizione per Ip costante, metto E0 costante
!         BT(LX,J) = BT(LX-1,J)*X2(LX-1)/X2(LX)
!     .        + E0/(ETAA(LX)*LX*X2(LX))
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
       BT(LX,J) = BT(LX-1,J)*
     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
     .        + E0*(2./LX)/(2.*X2(LX)*ETAA(LX)-VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
! nel caso J=0 non cambia nulla!

C     ...................B.C. PER BZ ALLA PARETE:                       
!         BZ(LX,J) = BZ(LX-1,J)
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
       BZ(LX,J) = BZ(LX-1,J)*
     .        (2.*ETAA(LX)+VRLX/LX)/
     .        (2.*ETAA(LX)-VRLX/LX)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
! nel caso J=0 non cambia nulla!

! daniele, 11/05/2007
! condizioni al bordo per vr finito
!         IF (IPR.EQ.-2) THEN
!            BZ(LX,J) = (2.*ETA0+ZDR*VR0(LX))/(2.*ETA0-ZDR*VR0(LX))*
!     .             BZ(LX-1,J)
!         END IF

C        BZ(LX,J) = BZ0(LX)                                             
         if (ippcd.eq.1) bzn(lx,j) =
     .          bzaini - it * ( bzaini-bzafin )/itend
         if (ippcd.eq.2) bzn(lx,j) =
     .          bzaini - ( bzaini-bzafin )*sin(4.*pi*it/itend)
         if (ippcd.eq.3) bzn(lx,j) = bzaini
         if (ippcd.eq.4) bzn(lx,j) =
     .          bzaini - ( bzaini-bzafin )*
     .          (1.-cos(4.*pi*it/itend))/2.

!       WRITE(6,*) 'bc2   0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0))) = ',
!     .        0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0)))

C     .............                                                     
72       IF(M.NE.0) GO TO 74                                            
         BZ(-1,J) = BZ(0,J)                                             
C  BR AUS BT +BZ                                                        
74     BR(0,J) =(0.0,0.0)                                               
       DO 150 IX=1,LX-1                                                 
       BR(IX,J)=(I/ZS22(IX-1)*(-BT(IX-1,J)*M-BZ(IX-1,J)*N*X2(IX-1)      
     1 /RR)+BR(IX-1,J)*X(IX-1))/X(IX)                                   
  150 CONTINUE                                                          
      IF(M.NE.1) GO TO 73                                               
C      BR(0,J) = 0.5*(BR(1,J)-BT(1,J))                                  
       BR(0,J) = (4.0*BR(1,J)-BR(2,J))*0.333333333                      
       BT(-1,J) = -BT(0,J)+2.0*I*BR(0,J)                                
   73 IF(J.EQ.0) GO TO 82                                               
       DO 40 IX=0,LX                                                    
       BR(IX,J1) = CONJG(BR(IX,J))                                      
   40 CONTINUE                                                          
       DO 402 IX=-1,LX                                                  
       BT(IX,J1) = CONJG(BT(IX,J))                                      
       BZ(IX,J1) = CONJG(BZ(IX,J))                                      
  402 CONTINUE                                                          
   82 CONTINUE                                                          
   28 CONTINUE                                                          

!       WRITE(6,*) 'dopo28   0.5*(REAL(BT(LX-1,0))+REAL(BT(LX,0))) = ',
!     .        0.5*(REAL(BT(LX-1,0))+REAL(BT(LX,0)))
                                                                        
       DO 41 IM=0,MY                                                    
       M = MM(IM)                                                       
       DO 41 N=NA1(M),NE(M)                                             
       J = JANZ(M,N)                                                    
       J1= JANZ(-M,-N)                                                  
         DO 43 IX=0,LX-1                                                
              R3(IX) = BZ(IX,J) - DT*ZDR*(ETAA(IX+1)-ETAA(IX))*0.5*       
     &                 I*N/RR*(BR(IX,J) + BR(IX+1,J))                   
 43      CONTINUE                                                       
                                                                        
         DO 901 L=1,2                                                   
              FF(L,0) = (0.0,0.0)                                       
         DO 901 K=1,2                                                   
              EE(L,K,0) = (0.0,0.0)                                     
 901     CONTINUE                                                       
                                                                        
         A12 = ETA2(0)*DT                                               
         A13 = DT*ZDR*0.5*(ETAA(1) - ETAA(0))                             
         KMN2 = (M/X2(0))**2 + (N/RR)**2                                
         ALFS = 1. + A12*ZDR**2*X(1)/X2(0) + A12*(KMN2+1./X2(0)**2)     
     &             - A13*ZDR*X(1)/X2(0)                                 
         BETS = 1. + A12*ZDR**2*X(1)/X2(0) + A12*KMN2 - A13*ZDR         
         IF (M.EQ.0)   BETS = BETS + 2.*A13*ZDR                         
         EE(2,1,0) = I*M/ALFS/X2(0)*(A12/X2(0) - A13)                   
         EE(2,2,0) = ZDR*X(1)/ALFS/X2(0)*(A12*ZDR + A13)                
         IF (M.EQ.1) THEN                                               
              EE(1,1,0) = (1.0,0.0)                                     
              EE(2,1,0) = 2.*EE(2,1,0)                                  
         END IF                                                         
         FF(2,0) = BT(0,J)/ALFS                                         
         AAA(0) = ZDR/BETS*(A12*ZDR*X(1)/X2(0) + A13)                   
         BBB(0) = R3(0)/BETS                                            
                                                                        
         DO 902 IX=1,LX-1                                               
              KMN1 = (M/X(IX))**2 + (N/RR)**2                           
              KMN2 = (M/X2(IX))**2 + (N/RR)**2                          
              A11 = ETAA(IX)*DT                                          
              A12 = ETA2(IX)*DT                                         
              A13 = DT*ZDR*0.5*(ETAA(IX+1) - ETAA(IX))                    
                                                                        
           AA(1,1,IX) = -A11*(ZDR**2 + ZDR*0.5/X(IX))                   
           AA(1,2,IX) = (0.0,0.0)                                       
           AA(2,1,IX) = -I*M/X2(IX)*(A12/X2(IX) - A13)                  
           AA(2,2,IX) = -ZDR*X(IX+1)/X2(IX)*(A12*ZDR + A13)             
                                                                        
           BB(1,1,IX) = 1. + A11*(2.*ZDR**2 + (KMN1 + 1./X(IX)**2))     
           BB(1,2,IX) =  A11*I*M/X(IX)**2                               
           BB(2,1,IX) = -I*M/X2(IX)*(A12/X2(IX) - A13)                  
           BB(2,2,IX) = 1. + A12*(2.*ZDR**2 + (KMN2 + 1./X2(IX)**2))    
     &                     - A13/X2(IX)                                 
                                                                        
           CC(1,1,IX) = -A11*(ZDR**2 - ZDR*0.5/X(IX))                   
           CC(1,2,IX) =  A11*I*M/X(IX)**2                               
           CC(2,1,IX) = (0.0,0.0)                                       
           CC(2,2,IX) = -ZDR*X(IX)/X2(IX)*(A12*ZDR - A13)               
                                                                        
           ALJ(IX)  = -A12*ZDR**2*X(IX+1)/X2(IX) - A13*ZDR              
           BETJ(IX) = 1. + A12*2.*ZDR**2 + A12*KMN2                     
           GAMJ(IX) = -A12*ZDR**2*X(IX)/X2(IX) + A13*ZDR                
                                                                        
           DD(1,IX) = BR(IX,J)                                          
           DD(2,IX) = BT(IX,J)                                          
                                                                        
              DO 903 L=1,2                                              
                   AV1(L) = FF(L,IX-1)                                  
                   AV2(L) = DD(L,IX)                                    
              DO 903 K=1,2                                              
                   AM1(L,K) = AA(L,K,IX)                                
                   AM2(L,K) = BB(L,K,IX)                                
                   AM3(L,K) = CC(L,K,IX)                                
                   AM4(L,K) = EE(L,K,IX-1)                              
 903          CONTINUE                                                  
              CALL MATCAL(AM1,AM2,AM3,AM4,AM5,AV1,AV2,AV3)              
              DO 904 L=1,2                                              
                   FF(L,IX) = AV3(L)                                    
              DO 904 K=1,2                                              
                   EE(L,K,IX) = AM5(L,K)                                
 904          CONTINUE                                                  
              AAA(IX) = -ALJ(IX)/(BETJ(IX) + GAMJ(IX)*AAA(IX-1))        
              BBB(IX) = (R3(IX) - GAMJ(IX)*BBB(IX-1))/                  
     &                           (BETJ(IX) + GAMJ(IX)*AAA(IX-1))        
 902     CONTINUE                                                       
                                                                        
!         RRAD = X2(LX-1)/(2. - X2(LX-1))
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
         RRAD = (2.*ETAA(LX)*X2(LX-1) + VRLX/LX)/
     .          (2.*ETAA(LX)*X2(LX) - VRLX/LX)
! Daniele, 23/06/2009
! questo è per imporre il campo magnetico radiale al bordo=0.
         UU(1,LX) = (0.0,0.0)
! Daniele, 18/08/2009
! caso generale per BR =/ 0
! ATTENZIONE! per ora metto J=1 perchè penso a simulazioni 2d
!         IF (J.EQ.1) UU(1,LX) = (0.0,1.E-4)
!         IF (J.EQ.1) UU(1,LX) = (0.0,5.E-4)
!         IF (J.EQ.1) UU(1,LX) = (0.0,1.E-3)
!!         IF (J.EQ.1) UU(1,LX) = (0.0,2.5E-3)
!         IF (J.EQ.1) UU(1,LX) = (0.0,2.E-3)
!         IF (J.EQ.1) UU(1,LX) = (0.0,5.E-3)
!         IF (J.EQ.1) UU(1,LX) = (0.0,1.E-2)
!!         IF (J.EQ.1) UU(1,LX) = (0.0,7.5E-3)
!         IF (J.EQ.1) UU(1,LX) = (0.0,1.E-2)
!         IF (J.EQ.1) UU(1,LX) = (0.0,1.5E-2)
!!         IF (J.EQ.1) UU(1,LX) = (0.0,2.E-2)
!!         IF (J.EQ.1) UU(1,LX) = (0.0,4.5E-2)

! ATTENZIONE! questo invece nel caso 3d
         IF (J.EQ.72) UU(1,LX) = (0.0,1.5E-4)

         IF (J.EQ.0) THEN                                               
!              UU(2,LX) = BT0(LX)                                        
!              UU(1,LX-1) = EE(1,2,LX-1)*UU(2,LX) + FF(1,LX-1)           
!              UU(2,LX-1) = EE(2,2,LX-1)*UU(2,LX) + FF(2,LX-1)           
! Daniele, 25/11/2008
! più giusto di quello prima (uso BTWALL invece di BT0)
!              UU(2,LX-1) = (2.*BTWALL*EE(2,2,LX-1) + FF(2,LX-1))/
!     .          (1.+EE(2,2,LX-1))
!              UU(2,LX) = 2.*BTWALL-UU(2,LX-1)
!              UU(1,LX-1) = EE(1,2,LX-1)*UU(2,LX) + FF(1,LX-1)
! Daniele, 25/11/2008
! invece di mettere condizione per Ip costante, metto E0 costante
!              UU(2,LX-1) = (EE(2,2,LX-1)*E0/(ETAA(LX)*X2(LX)*LX)
!     .             + FF(2,LX-1))/(1. - EE(2,2,LX-1)*RRAD)          
!              UU(2,LX) = UU(2,LX-1)*RRAD + E0/(ETAA(LX)*X2(LX)*LX)
!              UU(1,LX-1) = EE(1,2,LX-1)*UU(2,LX) + FF(1,LX-1)
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
!              UU(2,LX-1) = ((EE(2,2,LX-1)*E0*(2./LX))
!     .             /(2.*ETAA(LX)*X2(LX) - VRLX/LX) + FF(2,LX-1))
!     .             /(1. - EE(2,2,LX-1)*RRAD)          
!              UU(2,LX) = UU(2,LX-1)*RRAD +
!     .               E0*(2./LX)/(2.*ETAA(LX)*X2(LX) - VRLX/LX)
!              UU(1,LX-1) = EE(1,2,LX-1)*UU(2,LX) + FF(1,LX-1)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
              UU(2,LX-1) = (EE(2,1,LX-1)*UU(1,LX)+
     .               (EE(2,2,LX-1)*E0*(2./LX))
     .             /(2.*ETAA(LX)*X2(LX) - VRLX/LX) + FF(2,LX-1))
     .             /(1. - EE(2,2,LX-1)*RRAD)          
              UU(2,LX) = UU(2,LX-1)*RRAD +
     .               E0*(2./LX)/(2.*ETAA(LX)*X2(LX) - VRLX/LX)
              UU(1,LX-1) = EE(1,2,LX-1)*UU(2,LX) + FF(1,LX-1)

C      ..................B.C. PER BZ ALLA PARETE:                       
!              BZN(LX,J) = BBB(LX-1)/(1. - AAA(LX-1))
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
              BZN(LX,J) = BBB(LX-1)
     .               /((2.*ETAA(LX)-VRLX/LX)/(2.*ETAA(LX)+VRLX/LX)
     .               - AAA(LX-1))
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
! nel caso J=0 non cambia nulla!

!              BZN(LX-1,J) = BZN(LX,J)
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
              BZN(LX-1,J) = BZN(LX,J)/
     .        ((2.*ETAA(LX)+VRLX/LX)/
     .        (2.*ETAA(LX)-VRLX/LX))
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
! nel caso J=0 non cambia nulla!

! daniele, 11/05/2007
! condizioni al bordo per vr finito
!              IF (IPR.EQ.-2) THEN
!                 BZN(LX-1,J) =
!     .                  (2.*ETA0-ZDR*VR0(LX))/(2.*ETA0+ZDR*VR0(LX))*
!     .                  BZN(LX,J)
!              END IF

              if (ippcd.eq.1) then
                 bzn(lx,j) = bzaini - it * ( bzaini-bzafin )/itend
                 BZN(LX-1,J) = AAA(LX-1)*BZN(LX,J) + BBB(LX-1)
              endif
              if (ippcd.eq.2) then
                 bzn(lx,j) =
     .                  bzaini - ( bzaini-bzafin )*sin(4.*pi*it/itend)
                 BZN(LX-1,J) = AAA(LX-1)*BZN(LX,J) + BBB(LX-1)
              endif
              if (ippcd.eq.3) then
                 bzn(lx,j) = bzaini
                 BZN(LX-1,J) = AAA(LX-1)*BZN(LX,J) + BBB(LX-1)
              endif
              if (ippcd.eq.4) then
                 bzn(lx,j) =
     .                  bzaini - ( bzaini-bzafin )*
     .                  (1.-cos(4.*pi*it/itend))/2.
                 BZN(LX-1,J) = AAA(LX-1)*BZN(LX,J) + BBB(LX-1)
              endif


C      .................                                                
         ELSE                                                           
!              UU(2,LX-1) = FF(2,LX-1)/(1. - EE(2,2,LX-1)*RRAD)          
!              UU(1,LX-1) = EE(1,2,LX-1)*RRAD*UU(2,LX-1) + FF(1,LX-1)    
!              UU(2,LX) = UU(2,LX-1)*RRAD                                
! Daniele, 28/11/2008
! metto le B.C.s nuove anche per i modi diversi da zero
! mantenendo però nulle le componenti (m,n)=/(0,0) di vr(a)
! (questa è una differenza rispetto a PIXIE3D)
! già cambiare la def di RRAD sopra sistema quasi tutto
!              UU(2,LX-1) = FF(2,LX-1)/(1. - EE(2,2,LX-1)*RRAD)          
!              UU(2,LX) = UU(2,LX-1)*RRAD
!              UU(1,LX-1) = EE(1,2,LX-1)*UU(2,LX) + FF(1,LX-1)
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
              UU(2,LX-1) = (EE(2,1,LX-1)*UU(1,LX)+
     .             EE(2,2,LX-1)*(I*ETAA(LX)*M*UU(1,LX))*(2./LX)
     .             /(2.*ETAA(LX)*X2(LX) - VRLX/LX) + FF(2,LX-1))
     .             /(1. - EE(2,2,LX-1)*RRAD)          
              UU(2,LX) = UU(2,LX-1)*RRAD +
     .               I*ETAA(LX)*M*UU(1,LX)*(2./LX)
     .               /(2.*ETAA(LX)*X2(LX) - VRLX/LX)
              UU(1,LX-1) = EE(1,2,LX-1)*UU(2,LX) + FF(1,LX-1)

!              BZN(LX,J) = BBB(LX-1)/(1. - AAA(LX-1))                    
!              BZN(LX-1,J) = BZN(LX,J)                                   
! Daniele, 28/11/2008
! metto le B.C.s nuove anche per i modi diversi da zero
! mantenendo però nulle le componenti (m,n)=/(0,0) di vr(a)
! (questa è una differenza rispetto a PIXIE3D)
!              BZN(LX,J) = BBB(LX-1)
!     .               /((2.*ETAA(LX)-VRLX/LX)/(2.*ETAA(LX)+VRLX/LX)
!     .               - AAA(LX-1))
!              BZN(LX-1,J) = BZN(LX,J)/
!     .        ((2.*ETAA(LX)+VRLX/LX)/
!     .        (2.*ETAA(LX)-VRLX/LX))
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
              BZN(LX,J) = (BBB(LX-1)+(2./LX)
     .               /(2.*ETAA(LX) + VRLX/LX)
     .               *(I*ETAA(LX)*N/RR*UU(1,LX)))
     .               /((2.*ETAA(LX)-VRLX/LX)/(2.*ETAA(LX)+VRLX/LX)
     .               - AAA(LX-1))
              BZN(LX-1,J) = BZN(LX,J)/
     .               ((2.*ETAA(LX)+VRLX/LX)/
     .               (2.*ETAA(LX)-VRLX/LX))
     .               -(2./LX)/(2.*ETAA(LX) + VRLX/LX)
     .               *(I*ETAA(LX)*N/RR*UU(1,LX))

         END IF                                                         
                                                                        
         DO 907 IX=LX-2,0,-1                                            
              UU(1,IX) = EE(1,1,IX)*UU(1,IX+1) + EE(1,2,IX)*UU(2,IX+1)  
     &                                         + FF(1,IX)               
              UU(2,IX) = EE(2,1,IX)*UU(1,IX+1) + EE(2,2,IX)*UU(2,IX+1)  
     &                                         + FF(2,IX)               
              BZN(IX,J) = AAA(IX)*BZN(IX+1,J) + BBB(IX)                 
 907     CONTINUE                                                       
         DO 908 IX=0,LX                                                 
              BRN(IX,J) = UU(1,IX)                                      
              BTN(IX,J) = UU(2,IX)                                      
 908     CONTINUE                                                       
! Daniele, 10/05/2007
! per avere la condizione al bordo più giusta
! tenendo quella standard, non c'è bisogno di modifiche
!         IF (J.EQ.0) THEN
!            BTN(LX,J) = 2.*BTWALL-BTN(LX-1,J)
!            WRITE(6,*)
!     .             'dopo908: 0.5*(REAL(BTN(LX-1,0))+REAL(BTN(LX,0)))=',
!     .             0.5*(REAL(BTN(LX-1,0))+REAL(BTN(LX,0)))
!            WRITE(6,*)
!     .             'dopo908: REAL(BTN(LX,0))=',
!     .             REAL(BTN(LX,0))
!         END IF

         BTN(-1,J) = - BTN(0,J)                                         
         BZN(-1,J) = - BZN(0,J)                                         
         IF (M.NE.0) GO TO 909                                          
         BZN(-1,J) =   BZN(0,J)                                         
 909     IF (M.NE.1) GO TO 910                                          
         BTN(-1,J) = - BTN(0,J) + 2.*I*BRN(0,J)                         
 910     CONTINUE                                                       
                                                                        
C  BR AUS BT +BZ                                                        
       BRN(0,J) =(0.0,0.0)                                              
       DO 151 IX=1,LX-1                                                 
       BRN(IX,J)=(I/ZS22(IX-1)*(-BTN(IX-1,J)*M-BZN(IX-1,J)*N*X2(IX-1)   
     1 /RR)+BRN(IX-1,J)*X(IX-1))/X(IX)                                  
  151 CONTINUE                                                          
      IF(M.NE.1) GO TO 75                                               
       BRN(0,J) = (4.0*BRN(1,J)-BRN(2,J))*0.33333333333                 
       BTN(-1,J) = -BTN(0,J)+2.0*I*BRN(0,J)                             
                                                                        
   75 IF(J.EQ.0) GO TO 83                                               
       DO 45 IX=0,LX                                                    
       BRN(IX,J1) = CONJG(BRN(IX,J))                                    
   45 CONTINUE                                                          
       DO 404 IX=-1,LX                                                  
       BTN(IX,J1) = CONJG(BTN(IX,J))                                    
       BZN(IX,J1) = CONJG(BZN(IX,J))                                    
  404 CONTINUE                                                          
   83 CONTINUE                                                          
   41 CONTINUE                                                          
                                                                        
C..................  J-BESTIMMUNG ..................................... 
                                                                        
       DO 425 IM=0,MY                                                   
       M = MM(IM)                                                       
       DO 425 N=NA1(M),NE(M)                                            
       FAK1 = N/RR                                                      
       J = JANZ(M,N)                                                    
       J1 = JANZ(-M,-N)                                                 
       DO 426 IX=0,LX-1                                                 
       JR(IX,J) = I*(M/X2(IX)*BZN(IX,J)-FAK1*BTN(IX,J))                 
  426 CONTINUE                                                          
       DO 427 IX=1,LX-1                                                 
       JT(IX,J) = I*FAK1*BRN(IX,J)-2.0*ZS2(IX)*(BZN(IX,J)-BZN(IX-1,J))  
       JZ(IX,J) = 2.0*ZS2(IX)/X(IX)*(X2(IX)*BTN(IX,J)-X2(IX-1)          
     1 *BTN(IX-1,J))-I*M/X(IX)*BRN(IX,J)                                
  427 CONTINUE                                                          
       JT(LX,J) = (0.0,0.0)                                             
       JZ(LX,J) = (0.0,0.0)                                             
       JT(0,J) =(0.0,0.0)                                               
       JZ(0,J) =(0.0,0.0)                                               
       IF(M.NE.0) GO TO 422                                             
       JT(LX,J) = -2.*ZS2(LX)*(BZN(LX,J) - BZN(LX-1,J))                 
       JZ(LX,J) =  2.*ZS2(LX)*((2. - X2(LX-1))*BTN(LX,J) -              
     &                                          X2(LX-1)*BTN(LX-1,J))   
        JZ(0,J) = 2.*BTN(0,J)/X2(0)                                     
  422 IF(M.NE.1) GO TO 440                                              
       JT(0,J) = JT(1,J)                                                
  440 IF(J.EQ.0) GO TO 425                                              
       DO 439 IX=0,LX                                                   
       JR(IX,J1) = CONJG(JR(IX,J))                                      
  439 CONTINUE                                                          
       DO 419 IX=0,LX                                                   
       JT(IX,J1) = CONJG(JT(IX,J))                                      
       JZ(IX,J1) = CONJG(JZ(IX,J))                                      
  419 CONTINUE                                                          
  425 CONTINUE                                                          
C                                                                       
                                                                        
C ............................DENSITY ..................................
C                                                                       
C     DO 720 IM=0,MY                                                    
C      M = MM(IM)                                                       
C      DO 720 N=NA1(M),NE(M)                                            
C                                                                       
C       DO 721 IX = 0,LX                                                
C       SUM0(IX) = (0.0,0.0)                                            
C       SUM1(IX) = (0.0,0.0)                                            
C721     CONTINUE                                                       
C                                                                       
C       DO 722 MS = -MM(MY)+M,MM(MY)                                    
C      IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) GO TO 722                 
C        DO 723 NS = NA(MS),NE(MS)                                      
C      IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) GO TO 723               
C        J = JANZ(MS,NS)                                                
C        J1= JANZ(M-MS,N-NS)                                            
C          DO 724 IX=1,LX-1                                             
C                                                                       
C          VTIXJ= ((VT(IX-1,J)+VTN(IX-1,J))/2. +                        
C    .             (VT(IX,J)  +VTN(IX,J)  )/2.  )/2.                    
C          VZIXJ= ((VZ(IX-1,J)+VZN(IX-1,J))/2. +                        
C    .             (VZ(IX,J)  +VZN(IX,J)  )/2.  )/2.                    
C                                                                       
C     DENSITY :                                                         
C          S0 = - ZDR * 0.5/ X(IX)                                      
C          S1 = X(IX+1)*DENN(IX+1,J1)*(VR(IX+1,J)+VRN(IX+1,J))/2.       
C          S2 = X(IX-1)*DENN(IX-1,J1)*(VR(IX-1,J)+VRN(IX-1,J))/2.       
C          S3 = DENN(IX,J1)                                             
C          S4 = I * M * VTIXJ  / X(IX)                                  
C          S5 = I * N * VZIXJ   / RR                                    
C          SUM0(IX) = S0 *(S1-S2) - S3 *(S4+S5) + SUM0(IX)              
C                                                                       
C    PRESSURE :                                                         
C          S0 = - ZDR * 0.5 / X(IX)                                     
C          S1 = X(IX+1)*PRNN(IX+1,J1)*(VR(IX+1,J)+VRN(IX+1,J))/2.       
C          S2 = X(IX-1)*PRNN(IX-1,J1)*(VR(IX-1,J)+VRN(IX-1,J))/2.       
C          S3 = PRNN(IX,J1)                                             
C          S4 = I * M * VTIXJ  / X(IX)                                  
C          S5 = I * N * VZIXJ   / RR                                    
C          S6 = -2./3. * PRNN(IX,J1) * (                                
C    .          1./X(IX) * (                                            
C    .     X(IX+1)* (VR(IX+1,J)+VRN(IX+1,J))/2.  -                      
C    .     X(IX-1)* (VR(IX-1,J)+VRN(IX-1,J))/2.    ) *0.5 *ZDR          
C    .       + I * MS /X(IX) * VTIXJ + I * NS /RR * VZIXJ    )          
C          SUM1(IX) = S0 *(S1-S2) - S3 *(S4+S5) + S6 + SUM1(IX)         
C                                                                       
C724       CONTINUE                                                     
C723       CONTINUE                                                     
C722       CONTINUE                                                     
C                                                                       
C          J = JANZ(MM(IM),N)                                           
C          J1= JANZ(-MM(IM),-N)                                         
C          DO 725 IX=1,LX-1                                             
C                                                                       
C          AMX = 0.2                                                    
C          DEN(IX,J)=AMX*(DE(IX+1,J)+DE(IX-1,J))/2.+ (1.-AMX)*DE(IX,J)  
C    .                +  DT * SUM0(IX)                                  
C                                                                       
C          PRN(IX,J)=AMX*(PR(IX+1,J)+PR(IX-1,J))/2.+ (1.-AMX)*PR(IX,J)  
C    .                +  DT * SUM1(IX)                                  
C                                                                       
C          DEN(IX,J1)= CONJG(DEN(IX,J))                                 
C          PRN(IX,J1)= CONJG(PRN(IX,J))                                 
C725       CONTINUE                                                     
C          IF ( MM(IM) .EQ. 0 ) THEN                                    
C          DEN(0,J)  = DEN(1,J)                                         
C          DEN(0,J1) = CONJG(DEN(0,J))                                  
C          PRN(0,J)  = PRN (1,J)                                        
C          PRN(0,J1) = CONJG(PRN(0,J))                                  
C          ELSE                                                         
C          DEN(0,J)  = 0.0                                              
C          DEN(0,J1) = CONJG(DEN(0,J))                                  
C          PRN(0,J)  = 0.0                                              
C          PRN(0,J1) = CONJG(PRN(0,J))                                  
C          ENDIF                                                        
C          IF( J .EQ. 0 )THEN                                           
C          DEN(LX,J) = DEN(LX-1,J)                                      
C          ELSE                                                         
C          DEN(LX,J) = 0.0                                              
C          ENDIF                                                        
C          DEN(LX,J1)=CONJG(DEN(LX,J))                                  
C         PRN(LX,J) = PRN(LX-1,J)                                       
C         PRN(LX,J) = PR(LX,J)                                          
C         PRN(LX,J1)=CONJG(PRN(LX,J))                                   
C720       CONTINUE                                                     
C                                                                       
C......................... correzione del flusso toroidale:
       FTOR = 1.

       if (ippcd.ne.0) go to 992

       IF(MOD(IT,ITOR) .EQ. 1 .AND. IT .GT. IAUS) THEN
c        IF( TORFLX .LT. TORF1*(1.-0.01) ) THEN
        IF( TORFLX .LT. TORF1*(1.d0-0.005d0) ) THEN

        TORF = 0.0
        DO 990 IX=0,LX-1
        TORF=TORF+ 2.* X2(IX)*REAL(BZ(IX,J0))/ZDR
990     CONTINUE
        FTOR = TORF1 / TORF
        WRITE(6,*) '  perso 0.5% del flusso rispetto a it= 0 '
        WRITE(6,*) '  correzione, it = ',it
        WRITE(6,*) '  correzione, FTOR = ',FTOR
        call flush(6)

c        IF(TORFLX .LE. torf1 ) THEN
c        ITEND = (ITP+1) * IAUS
c        WRITE(6,*) ' al prossimo giro mi fermo (it= ',IT
c        ENDIF

        ENDIF

! daniele, 11/05/2007
        IF( TORFLX .GT. TORF1*(1.d0+0.005d0) ) THEN

        TORF = 0.0
        DO 989 IX=0,LX-1
        TORF=TORF+ 2.* X2(IX)*REAL(BZ(IX,J0))/ZDR
 989  CONTINUE
        FTOR = TORF1 / TORF
        WRITE(6,*) '  guadagnato 0.5% del flusso rispetto a it= 0 '
        WRITE(6,*) '  correzione, it = ',it
        call flush(6)

c        IF(TORFLX .LE. torf1 ) THEN
c        ITEND = (ITP+1) * IAUS
c        WRITE(6,*) ' al prossimo giro mi fermo (it= ',IT
c        ENDIF

        ENDIF
       ENDIF


C......................... REINIZIALIZZAZIONE :                         
 992   DO 77 IM=0,MY                                                    
       M=MM(IM)                                                         
       DO 77 N=NA1(M),NE(M)                                             
       J=JANZ(M,N)                                                      
       J1=JANZ(-M,-N)                                                   
!       write(6,*) J,J1
       DO 403 IX=0,LX                                                   
          VR(IX,J)  = VRN(IX,J)                                         
          VR(IX,J1) = VRN(IX,J1)                                        
          BR(IX,J)  = BRN(IX,J)                                         
          BR(IX,J1) = BRN(IX,J1)                                        
  403  CONTINUE                                                         
       DO 405 IX=-1,LX                                                  
          VT(IX,J)  = VTN(IX,J)                                         
          VT(IX,J1) = VTN(IX,J1)                                        
          VZ(IX,J)  = VZN(IX,J)                                         
          VZ(IX,J1) = VZN(IX,J1)                                        
          BT(IX,J)  = BTN(IX,J)                                         
          BT(IX,J1) = BTN(IX,J1)                                        
          BZ(IX,J)  = BZN(IX,J)                                         
          BZ(IX,J1) = BZN(IX,J1)                                        
          IF(J .EQ.J0) BZ(IX,J)  = FTOR * BZN(IX,J)                     
          IF(J1.EQ.J0) BZ(IX,J1) = FTOR * BZN(IX,J1)                    
C         DE(IX,J)  = DEN(IX,J)                                         
C         DE(IX,J1) = DEN(IX,J1)                                        
C         PR(IX,J)  = PRN(IX,J)                                         
C         PR(IX,J1) = PRN(IX,J1)                                        
  405  CONTINUE                                                         
   77  CONTINUE                                                         

! Daniele, 28/11/2008
! faccio il calcolo completo, ma sempre nel caso con solo vr00 finita
       B2=(0.5*(BZ(LX,0)+BZ(LX-1,0))*0.5*(BZ(LX,0)+BZ(LX-1,0))
     .        + 0.5*(BT(LX,0)+BT(LX-1,0))*0.5*(BT(LX,0)+BT(LX-1,0)))
       DO J=1,IZ
          B2 = B2 + 2.*(
     .          0.5*(BZ(LX,J)+BZ(LX-1,J))*0.5*(BZ(LX,J)+BZ(LX-1,J))
     .        + 0.5*(BT(LX,J)+BT(LX-1,J))*0.5*(BT(LX,J)+BT(LX-1,J)))
       ENDDO

! Daniele, 25/11/2008
! imposto ogni volta VRLX per seguire l'impostazione su E0
       VRLX = -E0*0.5*(BT(LX,0)+BT(LX-1,0))/B2
! Daniele, 29/11/2008
! pezza se voglio tornare al caso vecchio con vr0(LX)=0
!       VRLX = 0.

C                                                                       
C.................. END OF MAIN CALCULATION ..........................  
C                                                                       
C......................................... GROWTH RATE                  
! daniele per ppcd, 05/09/2006
! ho spostato il calcolo del flusso toroidale in modo che lo faccia
! ad ogni time step. Serve per fare il calcolo di quanto abbassare bza
       TORFLX = 0.0                                                     
       DO 993 IX=0,LX-1                                                 
         TORFLX = TORFLX + 2.*X2(IX)*REAL(BZ(IX,J0))/ZDR                
 993   CONTINUE                                                         
       YFITP = ( REAL(BZ(LX-1,J0))+REAL(BZ(LX,J0)) )/(2*TORFLX)       
       YTH = ( REAL(BT(LX-1,J0))+REAL(BT(LX,J0)) )/(2*TORFLX)       
       YQITP = X2(0)/RR*REAL(BZ(0,J0))/REAL(BT(0,J0))                 

       IF(MOD(IT,IAUS).NE.0) GO TO 91                                   
       ITP = ITP+1                                                      
       XT(ITP) = ZEIT                                                   
       J = JANZ(MPLO(2),NPLO(2))                                        
       YT(ITP)=(REAL(VR(IG,J))**2-REAL(VRALT)**2+AIMAG(VR(IG,J))**2     
     1 -AIMAG(VRALT)**2)/(DT*(REAL(VR(IG,J))**2                         
     2 +REAL(VRALT)**2+AIMAG(VR(IG,J))**2+AIMAG(VRALT)**2))             
       YZ(ITP)=(REAL(VZ(IG,J))**2-REAL(VZALT)**2+AIMAG(VZ(IG,J))**2     
     1 -AIMAG(VZALT)**2)/(DT*(REAL(VZ(IG,J))**2                         
     2 +REAL(VZALT)**2+AIMAG(VZ(IG,J))**2+AIMAG(VZALT)**2))             
!       TORFLX = 0.0                                                     
!       DO 993 IX=0,LX-1                                                 
!         TORFLX = TORFLX + 2.*X2(IX)*REAL(BZ(IX,J0))/ZDR                
! 993   CONTINUE                                                         
!       YF(ITP) = ( REAL(BZ(LX-1,J0))+REAL(BZ(LX,J0)) )/(2*TORFLX)       
!       YTH     = ( REAL(BT(LX-1,J0))+REAL(BT(LX,J0)) )/(2*TORFLX)       
!       YQ(ITP) = X2(0)/RR*REAL(BZ(0,J0))/REAL(BT(0,J0))   
       YF(ITP) = YFITP
       YQ(ITP) = YQITP
      WRITE(6,*) '   REAL(BZ(LX-1,J0)) =',REAL(BZ(LX-1,J0)),
     .        '   REAL(BZ(LX,J0)) = ', REAL(BZ(LX,J0))
       WRITE(6,*) '   REAL(BT(LX-1,J0)) =',REAL(BT(LX-1,J0)),
     .        '   REAL(BT(LX,J0)) = ', REAL(BT(LX,J0))
!       WRITE(6,*) '   0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0))) = ',
!     .        0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0)))
      WRITE(6,*) '   TH =',YTH, '   F = ', YF(ITP),'    ITP=', ITP      
      WRITE(6,*) '    flusso toroidale= ', TORFLX,'  ftor=',FTOR        
      WRITE(6,*) ' Qo= ', YQ(ITP),'    ITP=', ITP
      call flush(6)
C -  MAGNETIC AND KINETIC ENERGY                                        
C         Y0BR(ITP) = 0.0                                               
C         Y0VR(ITP) = 0.0                                               
C         DO 997 IX=0,LX-1                                              
C              Y0BR1 = (ABS(AIMAG(BR(IX,J))))**2 + (ABS(REAL(BT(IX,J)   
C     &                )))**2 + (ABS(REAL(BZ(IX,J))))**2                
C              Y0BR2 = (ABS(AIMAG(BR(IX+1,J))))**2 + (ABS(REAL(BT(IX+1,J
C     &                )))**2 + (ABS(REAL(BZ(IX+1,J))))**2              
C              Y0VR1 = (ABS(REAL(VR(IX,J))))**2 + (ABS(AIMAG(VT(IX,J)   
C     &                )))**2 + (ABS(AIMAG(VZ(IX,J))))**2               
C              Y0VR2 = (ABS(REAL(VR(IX+1,J))))**2 + (ABS(AIMAG(VT(IX+1,J
C     &                )))**2 + (ABS(AIMAG(VZ(IX+1,J))))**2             
C              Y0BR(ITP) = Y0BR(ITP) + 0.5*(Y0BR1 + Y0BR2)*X2(IX)/ZDR   
C              Y0VR(ITP) = Y0VR(ITP) + 0.5*(Y0VR1 + Y0VR2)*X2(IX)/ZDR   
C 997     CONTINUE                                                      
C         Y0BR(ITP) = 4.*PI**2*RR*Y0BR(ITP)                             
C         Y0VR(ITP) = 4.*PI**2*RR*Y0VR(ITP)                             
C         IF (MY.EQ.1) THEN                                             
C              Y1BR(ITP) = Y0BR(ITP)                                    
C              Y1VR(ITP) = Y0VR(ITP)                                    
C              GO TO 91                                                 
C         END IF                                                        
C         Y1BR(ITP) = 0.0                                               
C         Y1VR(ITP) = 0.0                                               
C         DO 995 IM=2,MY                                                
C              M = MM(IM)                                               
C              N = NANF(IM)                                             
C              J = JANZ(M,N)                                            
C         DO 995 IX=1,LX-1                                              
C              Y0BR1 = (ABS(AIMAG(BR(IX,J))))**2 + (ABS(REAL(BT(IX,J)   
C     &                )))**2 + (ABS(REAL(BZ(IX,J))))**2                
C              Y0BR2 = (ABS(AIMAG(BR(IX+1,J))))**2 + (ABS(REAL(BT(IX+1,J
C     &                )))**2 + (ABS(REAL(BZ(IX+1,J))))**2              
C              Y0VR1 = (ABS(REAL(VR(IX,J))))**2 + (ABS(AIMAG(VT(IX,J)   
C     &                )))**2 + (ABS(AIMAG(VZ(IX,J))))**2               
C              Y0VR2 = (ABS(REAL(VR(IX+1,J))))**2 + (ABS(AIMAG(VT(IX+1,J
C     &                )))**2 + (ABS(AIMAG(VZ(IX+1,J))))**2             
C              Y1BR(ITP) = Y1BR(ITP) + 0.5*(Y0BR1 + Y0BR2)*X2(IX)/ZDR   
C              Y1VR(ITP) = Y1VR(ITP) + 0.5*(Y0VR1 + Y0VR2)*X2(IX)/ZDR   
C 995     CONTINUE                                                      
C         Y1BR(ITP) = 4.*PI**2*RR*Y1BR(ITP) + Y0BR(ITP)                 
C         Y1VR(ITP) = 4.*PI**2*RR*Y1VR(ITP) + Y0VR(ITP)                 
C....................................... ENERGY CONTENTS                
                                                                        
        YBRITP=0.0
        YVRITP=0.0
        YBRFITP=0.0
        YVRFITP=0.0
        DO 995 IM = 0,MY                                                
         M = MM(IM)                                                     
         DO 995 N = NA1(M) , NE(M)                                      
          J = JANZ(M,N)                                                 
                                                                        
C          IF (ITP .EQ. 1) WRITE(6,*)                                   
C    .            '   M,N,J,ZDR',M,'  ',N,'  ',JANZ(M,N),'  ',ZDR       
           DO 995 IX=1,LX-1                                             
            YBR1= ((ABS(AIMAG(BR(IX,J))))**2. +                         
     .            (ABS(REAL (BT(IX,J))))**2. +                          
     .            (ABS(REAL (BZ(IX,J))))**2. )/2.                       
            YBR2= ((ABS(AIMAG(BR(IX+1,J))))**2. +                       
     .            (ABS(REAL (BT(IX+1,J))))**2. +                        
     .            (ABS(REAL (BZ(IX+1,J))))**2. )/2.                     
            YVR1= ((ABS(REAL (VR(IX,J))))**2. +                         
     .            (ABS(AIMAG(VT(IX,J))))**2. +                          
     .            (ABS(AIMAG(VZ(IX,J))))**2. )/2.                       
            YVR2= ((ABS(REAL (VR(IX+1,J))))**2. +                       
     .            (ABS(AIMAG(VT(IX+1,J))))**2. +                        
     .            (ABS(AIMAG(VZ(IX+1,J))))**2. )/2.                     
                                                                        
            IF( M .EQ. 0 .AND. N .EQ. 0 ) THEN
               FSUSI = 1.
            ELSE
               FSUSI = 2.
               YBRFITP = YBRFITP + 0.5 *FSUSI* (YBR1+YBR2)*X2(IX)/ZDR
               YVRFITP = YVRFITP + 0.5 *FSUSI* (YVR1+YVR2)*X2(IX)/ZDR
            ENDIF
            YBRITP = YBRITP + 0.5 *FSUSI* (YBR1+YBR2)*X2(IX)/ZDR
            YVRITP = YVRITP + 0.5 *FSUSI* (YVR1+YVR2)*X2(IX)/ZDR
995     CONTINUE                                                        
                                                                        
        YBTOT = 4.*PI**2.*RR* YBRITP                                    
        YVTOT = 4.*PI**2.*RR* YVRITP                                    
        YBF = 4.*PI**2.*RR* YBRFITP
        YVF = 4.*PI**2.*RR* YVRFITP
       WRITE(6,*) ' energia totale mag:',YBTOT
       WRITE(6,*) ' energia fluttuazioni mag:',YBF
       WRITE(6,*) ' energia totale kin:',YVTOT
       WRITE(6,*) ' energia fluttuazioni kin:',YVF
       call flush(6)

C........................................ WRITE ON DISK                 
   91  IF(MOD(IT,IB).NE.0) GO TO 117                                    
       IBIB  = IBIB + 1                                                 
                                                                        
C        IF ( IT .EQ. ITEND )                                           
C    .   REWIND IBDE                                                    
C        IF( IT .EQ. ITEND )                                            
C    .   WRITE(IBDE )      ZEIT,LX,IZ                                   
C        IF ( IT .EQ. ITEND )                                           
C    .   REWIND IBPR                                                    
C        IF( IT .EQ. ITEND )                                            
C    .   WRITE(IBPR )      ZEIT,LX,IZ                                   
                                                                        
         IF ( IT .EQ. ITEND )                                           
     .   REWIND IBAUS                                                   
         IF( IT .EQ. ITEND )                                            
     .   WRITE(IBAUS)      ZEIT,LX,IZ                                   
                                                                        
C      WRITE(IBDE1)   ZEIT,LX,IZ                                        
C      WRITE(IBPR1)   ZEIT,LX,IZ                                        
       WRITE(IBOUT1)  ZEIT,LX,IZ                                        
       DO 252 J=0,IZ                                                    
          AVT = AIMAG(VT(-1,J))                                         
          AVZ = AIMAG(VZ(-1,J))                                         
          ABT =  REAL(BT(-1,J))                                         
          ABZ =  REAL(BZ(-1,J))                                         
                                                                        
          IF ( IT .EQ. ITEND )                                          
     .    WRITE(IBAUS)      AVT, AVZ, ABT, ABZ                          
                                                                        
          WRITE(IBOUT1)     AVT, AVZ, ABT, ABZ                          
                                                                        
       DO 252 IX=0,LX                                                   
          AVR =  REAL(VR(IX,J))                                         
          AVT = AIMAG(VT(IX,J))                                         
          AVZ = AIMAG(VZ(IX,J))                                         
          ABR = AIMAG(BR(IX,J))                                         
          ABT =  REAL(BT(IX,J))                                         
          ABZ =  REAL(BZ(IX,J))                                         
C         ADE =  REAL(DE(IX,J))                                         
C         APR =  REAL(PR(IX,J))                                         
                                                                        
C         IF ( IT .EQ. ITEND ) WRITE(IBDE )   ADE                       
C         IF ( IT .EQ. ITEND ) WRITE(IBPR )   APR                       
C         WRITE(IBDE1)   ADE                                            
C         WRITE(IBPR1)   APR                                            
                                                                        
          IF (IT .EQ. ITEND) WRITE(IBAUS) AVR, AVT, AVZ, ABR, ABT, ABZ  
          WRITE(IBOUT1) AVR, AVT, AVZ, ABR, ABT, ABZ                    
                                                                        
 252   CONTINUE                                                         

       call flush(IBAUS)
       call flush(IBOUT1)
 
C      IF(IT .EQ. ITEND) THEN                                           
C                                                                       
C      WRITE(6,*)'  zeit ', ZEIT,LX,IZ                                  
C      WRITE(6,*)' vr  oo =',VR(20,JANZ(0,0)),'  ',VR(50,JANZ(0,0))     
C      WRITE(6,*)' vt  oo =',VT(20,JANZ(0,0)),'  ',VT(50,JANZ(0,0))     
C      WRITE(6,*)' vz  oo =',VZ(20,JANZ(0,0)),'  ',VZ(50,JANZ(0,0))     
C      WRITE(6,*)' br  oo =',BR(20,JANZ(0,0)),'  ',BR(50,JANZ(0,0))     
C      WRITE(6,*)' bt  oo =',BT(20,JANZ(0,0)),'  ',BT(50,JANZ(0,0))     
C      WRITE(6,*)' bz  oo =',BZ(20,JANZ(0,0)),'  ',BZ(50,JANZ(0,0))     
C      WRITE(6,*)'    janzoo',JANZ(0,0)                                 
C      WRITE(6,*)' vr  11 =',VR(20,JANZ(1,1)),'  ',VR(50,JANZ(1,1))     
C      WRITE(6,*)' vt  11 =',VT(20,JANZ(1,1)),'  ',VT(50,JANZ(1,1))     
C      WRITE(6,*)' vz  11 =',VZ(20,JANZ(1,1)),'  ',VZ(50,JANZ(1,1))     
C      WRITE(6,*)' br  11 =',BR(20,JANZ(1,1)),'  ',BR(50,JANZ(1,1))     
C      WRITE(6,*)' bt  11 =',BT(20,JANZ(1,1)),'  ',BT(50,JANZ(1,1))     
C      WRITE(6,*)' bz  11 =',BZ(20,JANZ(1,1)),'  ',BZ(50,JANZ(1,1))     
C      WRITE(6,*)' janz11,23  ',JANZ(1,1),'  ',JANZ(2,3)                
C                                                                       
C      WRITE(6,*)'   '                                                  
C      WRITE(6,*)' vr  23 =',VR(20,JANZ(2,3)),'  ',VR(50,JANZ(2,3))     
C      WRITE(6,*)' vt  23 =',VT(20,JANZ(2,3)),'  ',VT(50,JANZ(2,3))     
C      WRITE(6,*)' vz  23 =',VZ(20,JANZ(2,3)),'  ',VZ(50,JANZ(2,3))     
C      WRITE(6,*)' br  23 =',BR(20,JANZ(2,3)),'  ',BR(50,JANZ(2,3))     
C      WRITE(6,*)' bt  23 =',BT(20,JANZ(2,3)),'  ',BT(50,JANZ(2,3))     
C      WRITE(6,*)' bz  23 =',BZ(20,JANZ(2,3)),'  ',BZ(50,JANZ(2,3))     
C      WRITE(6,*)'   '                                                  
C      WRITE(6,*)' vr -1-1 =',VR(20,JANZ(-1,-1)),'  ',VR(50,JANZ(-1,-1))
C      WRITE(6,*)' vt -1-1 =',VT(20,JANZ(-1,-1)),'  ',VT(50,JANZ(-1,-1))
C      WRITE(6,*)' vz -1-1 =',VZ(20,JANZ(-1,-1)),'  ',VZ(50,JANZ(-1,-1))
C      WRITE(6,*)' br -1-1 =',BR(20,JANZ(-1,-1)),'  ',BR(50,JANZ(-1,-1))
C      WRITE(6,*)' bt -1-1 =',BT(20,JANZ(-1,-1)),'  ',BT(50,JANZ(-1,-1))
C      WRITE(6,*)' bz -1-1 =',BZ(20,JANZ(-1,-1)),'  ',BZ(50,JANZ(-1,-1))
C      WRITE(6,*) '  '                                                  
C         ENDIF                                                         
                                                                        
CC       WRITE(IBAUS) ITP,(YT(IX),IX=1,ITP),(YZ(IX),IX=1,ITP),          
CC     1   (YF(IX),IX=1,ITP),  (YQ(IX),IX=1,ITP),                       
CC     2 (Y0BR(IX),IX=1,ITP),(Y1BR(IX),IX=1,ITP),                       
CC     3 (Y0VR(IX),IX=1,ITP),(Y1VR(IX),IX=1,ITP),(XT(IX),IX=1,ITP)      
CC       REWIND IBOUT2                                                  
CC       WRITE(IBOUT2,101) ITP                                          
CC       WRITE(IBOUT2,102) (XT(IX),IX=1,ITP)                            
CC       WRITE(IBOUT2,102)   (YF(IX),IX=1,ITP),  (YQ(IX),IX=1,ITP)      
CC       WRITE(IBOUT2,102)   (YT(IX),IX=1,ITP),  (YZ(IX),IX=1,ITP)      
CC       WRITE(IBOUT2,102) (Y0BR(IX),IX=1,ITP),(Y1BR(IX),IX=1,ITP)      
CC       WRITE(IBOUT2,102) (Y0VR(IX),IX=1,ITP),(Y1VR(IX),IX=1,ITP)      
 117   CONTINUE                                                         
c       titr = second() - tt
c       ttot = ttot + titr
c       tt = second()
c       print*,'            cpu     ',titr,' sec '
c       print*,'            tot cpu ',ttot,' sec (',iz,'), it ',it
   50 CONTINUE                                                          

5000  continue

C................ MAIN LOOP ENDS HERE ................................  
c       t1 = second()
c       print*,'            MAIN LOOP ENDS HERE '
c       print*,'            tot cpu ',t1-t0,' (',ttot,
c     .                     ') sec (',iz,'), nt ',itend
c                                                                        
c1111     WRITE(6,1002) X2(LX-3),BZ(LX-3,0),BT(LX-3,0)                   
c         WRITE(6,1002) X2(LX-2),BZ(LX-2,0),BT(LX-2,0)                   
c         WRITE(6,1002) X2(LX-1),BZ(LX-1,0),BT(LX-1,0)                   
c         WRITE(6,1002) X2(LX),BZ(LX,0),BT(LX,0)                         
c         WRITE(6,550) ITEND,ITP,ZEIT,RS                                 
1111         write(6,*) ' zeit= ',zeit,'  stato finale:'
         WRITE(6,*) 'x= ',X2(LX-4),' bz ',BZ(LX-4,0),' bt ',BT(LX-4,0)
         WRITE(6,*) 'x= ',X2(LX-3),' bz ',BZ(LX-3,0),' bt ',BT(LX-3,0)
         WRITE(6,*) 'x= ',X2(LX-2),' bz ',BZ(LX-2,0),' bt ',BT(LX-2,0)
         WRITE(6,*) 'x= ',X2(LX-1),' bz ',BZ(LX-1,0),' bt ',BT(LX-1,0)
         WRITE(6,*) 'x= ',X2(LX),' bz ',BZ(LX,0),' bt ',BT(LX,0)
         call flush(6)

  100 FORMAT(6E12.4)                                                    
  101 FORMAT(I4)                                                        
  102 FORMAT(10E12.4)                                                   
  103 FORMAT(E12.4,2I4)                                                 
  550 FORMAT(1X,' ITEND,ITP=',2I10,' ZEIT=',1PE12.4,' RS=',E12.4)       
!       STOP                                                             
       END                                                              

      SUBROUTINE EQUIL (S,BZ,BT,P)
      DIMENSION Y0(3),DY0(3)
      COMMON /COMEQ/ Q01, AA1, BB1, RR, THETA0, ALPHA, CHI, BETA0, IPR
C
      IF (IPR.EQ.2) THETA0 = 1./(RR*Q01)
C     INITIAL CONDITIONS
       RADI = 1.E-8
      Y0(1) = 1.0
      Y0(2) = RADI*THETA0
      Y0(3) = BETA0
      NE = 3
      NP = 20
      RADF = S
      DRAD = (RADF - RADI)/(NP - 1)
      IF (S.EQ.0.) THEN
            BZ = 1.0
            BT = 0.0
             P = BETA0
            DBZ = 0.0
            DBT = THETA0
             DP = 0.0
            RETURN
      END IF
C     INTEGRATION
      DO 1 I=1,NP-1
            RAD0 = RADI + (I - 1)*DRAD
            CALL KMINT (NE,Y0,RAD0,DRAD,DY0)
 1    CONTINUE
      BZ = Y0(1)
      BT = Y0(2)
       P = Y0(3)
      DBZ = DY0(1)
      DBT = DY0(2)
       DP = DY0(3)
      RETURN
      END
C
      SUBROUTINE KMINT (NE,Y0,T0,DT,ER0)
      DIMENSION  Y0(3), YY(3), A(3), B(3), C(3), D(3), E(3), F(3),
     &          ER0(3)
      COMMON /COMEQ/ Q01, AA1, BB1, RR, THETA0, ALPHA, CHI, BETA0, IPR
C
      DO 1 K=1,NE
            YY(K) = Y0(K)
 1    CONTINUE
      T = T0
      DO 2 J=1,6
            RA = THETA0*RMU(T)
            RB = YY(1)*YY(1) + YY(2)*YY(2)
            RC = RA*RB/YY(2)
            RD = RC - YY(1)/T
            RE = CHI*T*RD*RD
            RF = RE/RB/2.
                  E(1) = - 2.*RA*YY(2) + RF*YY(1)
                  E(2) =   2.*RA*YY(1) + RF*YY(2) - YY(2)/T
                  E(3) = - RE
            GO TO (11,12,13,14,15,16), J
 11         DO 21 K=1,NE
            A(K) = DT*E(K)
 21         YY(K) = Y0(K) + A(K)/3.
            T = T0 + DT/3.
            GO TO 2
 12         DO 22 K=1,NE
            B(K) = DT*E(K)
 22         YY(K) = Y0(K) + A(K)/6. + B(K)/6.
            GO TO 2
 13         DO 23 K=1,NE
            C(K) = DT*E(K)
 23         YY(K) = Y0(K) + A(K)/8.+ C(K)*3./8.
            T = T0 + DT/2.
            GO TO 2
 14         DO 24 K=1,NE
            D(K) = DT*E(K)
 24         YY(K) = Y0(K) + A(K)/2. - C(K)*3./2. + 2.*D(K)
            T = T0 + DT
            GO TO 2
 15         DO 25 K=1,NE
            F(K) = DT*E(K)
 25         YY(K) = Y0(K) + A(K)/6. + D(K)*2./3. + F(K)/6.
            GO TO 2
 16         CONTINUE
 2    CONTINUE
      DO 3 K=1,NE
            Y0(K) = YY(K)
C           ER0(K) = (2.*A(K) - 9.*C(K) + 8.*D(K) - F(K))/3.
            ER0(K) = E(K)
 3    CONTINUE
      RETURN
      END
C
      FUNCTION RMU(X)
      COMMON /COMEQ/ Q01, AA1, BB1, RR, THETA0, ALPHA, CHI, BETA0, IPR
      IF (IPR.LT.2) THEN
c         RMU = 1. - IPR*X**ALPHA
      PI=3.141592653589793                                     
         RMU = ((1. + COS (PI*X))/2.)**ALPHA
! daniele, 08/04/2008
! non so perchè veniva usato il mu di prima, ora metto quello di mu&p
! aa1 viene usato come a nel caso voglio parete lontana dal plasma
!      IF (X.LE.AA1) THEN
!         RMU = 1./AA1 * (1.-(X/AA1)**ALPHA)
!      ELSE
!         RMU = 0.
!      ENDIF         
!         WRITE(6,*) 'x= ',X,' rmu ',RMU,' aa1= ',AA1,' alpha= ',alpha

      ELSE
         Q = Q01*(1. - AA1*X**2 + BB1*X**4)
         DQ = Q01*( - 2.*AA1*X + 4.*BB1*X**3)
         RMU = Q01/2.*(2.*Q - X*DQ)/((X/RR)**2 + Q**2)
      END IF
      RETURN
      END
 
         SUBROUTINE MATCAL(AM1,AM2,AM3,AM4,AM5,AV1,AV2,AV3)
         COMPLEX AM1(2,2), AM2(2,2), AM3(2,2), AM4(2,2), AM5(2,2),
     &      DET, AV1(2), AV2(2), AV3(2), BM1(2,2), BM2(2,2), BV1(2)
         BM1(1,1) = AM3(1,1)*AM4(1,1) + AM3(1,2)*AM4(2,1)
         BM1(1,2) = AM3(1,1)*AM4(1,2) + AM3(1,2)*AM4(2,2)
         BM1(2,1) = AM3(2,1)*AM4(1,1) + AM3(2,2)*AM4(2,1)
         BM1(2,2) = AM3(2,1)*AM4(1,2) + AM3(2,2)*AM4(2,2)
         BM1(1,1) = BM1(1,1) + AM2(1,1)
         BM1(1,2) = BM1(1,2) + AM2(1,2)
         BM1(2,1) = BM1(2,1) + AM2(2,1)
         BM1(2,2) = BM1(2,2) + AM2(2,2)
         DET = BM1(1,1)*BM1(2,2) - BM1(1,2)*BM1(2,1)
         BM2(1,1) =  BM1(2,2)/DET
         BM2(1,2) = -BM1(1,2)/DET
         BM2(2,1) = -BM1(2,1)/DET
         BM2(2,2) =  BM1(1,1)/DET
         BV1(1) = AV2(1) - AM3(1,1)*AV1(1) - AM3(1,2)*AV1(2)
         BV1(2) = AV2(2) - AM3(2,1)*AV1(1) - AM3(2,2)*AV1(2)
         AM5(1,1) = - BM2(1,1)*AM1(1,1) - BM2(1,2)*AM1(2,1)
         AM5(1,2) = - BM2(1,1)*AM1(1,2) - BM2(1,2)*AM1(2,2)
         AM5(2,1) = - BM2(2,1)*AM1(1,1) - BM2(2,2)*AM1(2,1)
         AM5(2,2) = - BM2(2,1)*AM1(1,2) - BM2(2,2)*AM1(2,2)
         AV3(1) = BM2(1,1)*BV1(1) + BM2(1,2)*BV1(2)
         AV3(2) = BM2(2,1)*BV1(1) + BM2(2,2)*BV1(2)
         RETURN
         END
 
       SUBROUTINE FCT(X,F,DERF)                                         
       COMMON /FCTX/CX1,CX2,ALC1,ALC2,TAX0,IXXX,NX1,CLX,RS              
       XX = X/CLX                                                       
       X1 = XX-RS                
C      F = NX1*(CX*TANH(ALC*X1 )+(1.0-CX)*XX+CX*TANH(ALC*RS))/TAX0-IXXX 
C      DERF = NX1/CLX*(CX*ALC/COSH(ALC *X1)**2+1.0-CX)/TAX0             
       F= NX1*(CX1*TANH(ALC1*XX)+CX2*(TANH(ALC2*X1)+TANH(ALC2*          
     1 RS))+(1.0-CX1-CX2)*XX)/TAX0-IXXX                                 
       DERF = NX1/CLX*(CX1*ALC1/COSH(ALC1*XX)**2                        
     1 +CX2*ALC2/COSH(ALC2*X1)**2+1.0-CX1-CX2)/TAX0                     
       RETURN                                                           
       END                                                              

      SUBROUTINE RTNI(X,F,DERF,FCT,XST,EPS,IEND,IER)                    
      IER=0                                                             
      X=XST                                                             
      TOL=X                                                             
      CALL FCT(TOL,F,DERF)                                              
      TOLF=100.*EPS                                                     
      DO 6 I=1,IEND                                                     
      IF(F)1,7,1                                                        
    1 IF(DERF)2,8,2                                                     
    2 DX=F/DERF                                                         
      X=X-DX                                                            
      TOL=X                                                             
      CALL FCT(TOL,F,DERF)                                              
      TOL=EPS                                                           
      A=ABS(X)                                                          
      IF(A-1.)4,4,3                                                     
    3 TOL=TOL*A                                                         
    4 IF(ABS(DX)-TOL)5,5,6                                              
    5 IF(ABS(F)-TOLF)7,7,6                                              
    6 CONTINUE                                                          
      IER=1                                                             
    7 RETURN                                                            
    8 IER=2                                                             
      RETURN                                                            
      END                                                               

      FUNCTION DISVAC2(x,dis0,alpha,rplasma)

      if (x.le.rplasma) then
         DISVAC2=dis0
      else
         DISVAC2=dis0*alpha
      end if

      END

      FUNCTION DISVAC3(x,dis0,alpha,beta,gamma,rplasma)

      DISVAC3=dis0*(1.+ (alpha-1.)*(min(x,rplasma)/rplasma)**beta*(
     .       (atan(gamma*(x-rplasma))-atan(-gamma*rplasma))/
     .       (atan(gamma*(1-rplasma))-atan(-gamma*rplasma))))
      END

      FUNCTION DISVAC4(x,dis0,alpha,beta,rplasma)

      if (x.le.rplasma) then
         DISVAC4=dis0*(1. + (alpha-1.)*(x/rplasma)**beta)
      else
         DISVAC4=dis0*alpha
      end if

      END

      include 'eq_csha.blc.for'
