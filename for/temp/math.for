      module math
       contains

!#! Marco, commento tutto, non mi serve!
!       SUBROUTINE EQUIL (S,BZ,BT,P)
!!#! Marco
!       implicit none
!       real(8) :: s,bz,bt,p
!       real(8) :: beta0,dbt,dbz,dp,drad,dy0
!       real(8) :: rad0,radf,radi,theta0,y0
!       integer :: i,ipr,ne,np
!
!       DIMENSION Y0(3),DY0(3)
!!      COMMON /COMEQ/ Q01, AA1, BB1, RR, THETA0, ALPHA, CHI, BETA0, IPR
!C
!!      IF (IPR.EQ.2) THETA0 = 1./(RR*Q01)
!C     INITIAL CONDITIONS
!       RADI = 1.E-8
!      Y0(1) = 1.0
!      Y0(2) = RADI*THETA0
!      Y0(3) = BETA0
!      NE = 3
!      NP = 20
!      RADF = S
!      DRAD = (RADF - RADI)/(NP - 1)
!      IF (S.EQ.0.) THEN
!            BZ = 1.0
!            BT = 0.0
!             P = BETA0
!            DBZ = 0.0
!            DBT = THETA0
!             DP = 0.0
!            RETURN
!      END IF
!C     INTEGRATION
!      DO 1 I=1,NP-1
!            RAD0 = RADI + (I - 1)*DRAD
!            CALL KMINT (NE,Y0,RAD0,DRAD,DY0)
! 1    CONTINUE
!      BZ = Y0(1)
!      BT = Y0(2)
!       P = Y0(3)
!      DBZ = DY0(1)
!      DBT = DY0(2)
!       DP = DY0(3)
!      RETURN
!      END
C

      SUBROUTINE KMINT (NE,Y0,T0,DT,ER0)
      DIMENSION  Y0(3), YY(3), A(3), B(3), C(3), D(3), E(3), F(3),
     &          ER0(3)
!      COMMON /COMEQ/ Q01, AA1, BB1, RR, THETA0, ALPHA, CHI, BETA0, IPR
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
!      COMMON /COMEQ/ Q01, AA1, BB1, RR, THETA0, ALPHA, CHI, BETA0, IPR
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
!         COMPLEX AM1(2,2), AM2(2,2), AM3(2,2), AM4(2,2), AM5(2,2),
!     &      DET, AV1(2), AV2(2), AV3(2), BM1(2,2), BM2(2,2), BV1(2)
!#! Marco, passo a complex(8)
         COMPLEX(8) AM1(2,2), AM2(2,2), AM3(2,2), AM4(2,2), AM5(2,2),
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
!#! Marco
      real(8) :: x

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

!       SUBROUTINE FCT(X,F,DERF)                                         
!       COMMON /FCTX/CX1,CX2,ALC1,ALC2,TAX0,IXXX,NX1,CLX,RS              
!       XX = X/CLX                                                       
!       X1 = XX-RS                                                       
!C      F = NX1*(CX*TANH(ALC*X1 )+(1.0-CX)*XX+CX*TANH(ALC*RS))/TAX0-IXXX 
!C      DERF = NX1/CLX*(CX*ALC/COSH(ALC *X1)**2+1.0-CX)/TAX0             
!       F= NX1*(CX1*TANH(ALC1*XX)+CX2*(TANH(ALC2*X1)+TANH(ALC2*          
!     1 RS))+(1.0-CX1-CX2)*XX)/TAX0-IXXX                                 
!       DERF = NX1/CLX*(CX1*ALC1/COSH(ALC1*XX)**2                        
!     1 +CX2*ALC2/COSH(ALC2*X1)**2+1.0-CX1-CX2)/TAX0                     
!       RETURN                                                           
!       END                                                              
!      SUBROUTINE RTNI(X,F,DERF,FCT,XST,EPS,IEND,IER)                    
!      IER=0                                                             
!      X=XST                                                             
!      TOL=X                                                             
!      CALL FCT(TOL,F,DERF)                                              
!      TOLF=100.*EPS                                                     
!      DO 6 I=1,IEND                                                     
!      IF(F)1,7,1                                                        
!    1 IF(DERF)2,8,2                                                     
!    2 DX=F/DERF                                                         
!      X=X-DX                                                            
!      TOL=X                                                             
!      CALL FCT(TOL,F,DERF)                                              
!      TOL=EPS                                                           
!      A=ABS(X)                                                          
!      IF(A-1.)4,4,3                                                     
!    3 TOL=TOL*A                                                         
!    4 IF(ABS(DX)-TOL)5,5,6                                              
!    5 IF(ABS(F)-TOLF)7,7,6                                              
!    6 CONTINUE                                                          
!      IER=1                                                             
!    7 RETURN                                                            
!    8 IER=2                                                             
!      RETURN                                                            
!      END                                                               
!#! Roba di Daniele
!       FUNCTION DISVAC2(x,dis0,alpha,rplasma)
!      
!       if (x.le.rplasma) then
!          DISVAC2=dis0
!       else
!          DISVAC2=dis0*alpha
!       end if
!      
!       END
!      
!       FUNCTION DISVAC3(x,dis0,alpha,beta,gamma,rplasma)
!      
!       DISVAC3=dis0*(1.+ (alpha-1.)*(min(x,rplasma)/rplasma)**beta*(
!     .       (atan(gamma*(x-rplasma))-atan(-gamma*rplasma))/
!     .       (atan(gamma*(1-rplasma))-atan(-gamma*rplasma))))
!       END
!      
!       FUNCTION DISVAC4(x,dis0,alpha,beta,rplasma)
!      
!       if (x.le.rplasma) then
!          DISVAC4=dis0*(1. + (alpha-1.)*(x/rplasma)**beta)
!       else
!          DISVAC4=dis0*alpha
!       end if
!
!       END
      end module math
