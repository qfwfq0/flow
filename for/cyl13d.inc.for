C  MY=numero degli M(senza NULL) 
C  N0=numero degli N (M=0) (senza NULL)
C  N1...N30= numero degli  N con tutti gli  MY (con il NULL)

       integer , parameter :: wp = kind(1.d0) ! working precision
       integer , parameter :: LX=100
       integer , parameter :: MY=2
!
!       scelta degli spettri:
!
       integer, parameter :: N0 = 25
       integer, parameter :: N1 = 65
       integer, parameter :: N2 = 45
       integer, parameter :: N3 = 0
       integer, parameter :: N4 = 0
       integer , parameter :: IZ=135
!       integer, parameter :: IZ = N0 + MY * N1
       integer , parameter :: numz = 300
       integer , parameter :: ntp= 4000

       real(wp) :: pi
       DATA   PI/3.14159265358979323846d0/
 
       INTEGER  :: MM(0:MY),NANF(0:MY),MWERT(-my:my)
       integer :: na1(0:my),na(0:my),ne(-my:my)
       integer :: nz(0:my),nzcon(0:my)
       integer :: n2d,m2d
       DATA MM/0,1,2/
       DATA NANF/-25,-55,-50/
!#! Marco, per non includere cyl2.inc.for
       DATA NZ/25,65,45/
       DATA NZCON/26,65,45/

!#! Marco, perturbazione iniziale con fasi sganciate
       logical ,parameter :: fasi_sganciate=.true.
!#! Marco, irrotational MP
       logical ,parameter :: irrot=.true.
       real(wp), parameter :: br_pert_freq=0.d0
       integer, parameter :: qmp = 0
!#! amplitude and phase
       real (wp) , dimension(0:qmp), parameter ::
     .             mp_perc=(/0.d0/)
       real (wp) , dimension(0:qmp), parameter :: 
     .             ph_mp=(/0.d0/)
!#! m_MP and n_MP (2D case)
       integer , dimension(0:qmp), parameter :: m_mp = (/0/)
       integer , dimension(0:qmp), parameter :: n_mp = (/0/)
       integer , dimension(0:qmp) :: j_mp, temp_mp
       real (wp) , dimension(0:qmp) :: phmp, amp_a
!#! old MP
       integer , dimension(0:qmp), parameter :: m_mpold = (/0/)
       integer , dimension(0:qmp), parameter :: n_mpold = (/0/)
       integer , dimension(0:qmp) :: j_mpold
       real (wp) , dimension(0:qmp), parameter :: 
     .             ampold=(/0.d0/)
       real (wp) :: pitol
       data   pitol/1.e-4/
    
!#! Marco, momentum source
       real(wp) :: m0_momsour(0:2)
       data m0_momsour/0.d0,0.d0,0.d0/
       complex(wp) ::momsour(0:lx,0:2)

!#! Marco, trick to break phase constanc
       COMPLEX(wp) , allocatable :: BR(:,:),BT(:,:)
     .             ,BZ(:,:)   
       COMPLEX(wp) , allocatable :: VR(:,:),VT(:,:)
     .             ,VZ(:,:)   
       COMPLEX(wp) , allocatable :: BRN(:,:),BTN(:,:)
     .             ,BZN(:,:)   
       COMPLEX(wp) , allocatable :: VRN(:,:),VTN(:,:)
     .             ,VZN(:,:) 
       COMPLEX(wp) , allocatable :: VRP(:,:),VTP(:,:)                      
       COMPLEX(wp) , allocatable :: JR(:,:),JT(:,:)
     .             ,JZ(:,:)   
       COMPLEX(wp) , allocatable :: ER(:,:), ET(:,:)
     .             ,EZ(:,:)
       COMPLEX(wp) , allocatable :: PHI(:,:),RHOC(:,:)  
! 27/11/03 Daniele: per il calcolo di Az e Atheta (nel gauge Ar=0)
	COMPLEX (wp) , allocatable :: AT(:,:),AZ(:,:)
! 27/02/04 Daniele: per il calcolo di Az e Atheta vecchi
	COMPLEX (wp) , allocatable :: ATOLD(:,:),AZOLD(:,:)

       REAL (wp) :: ZS2(0:LX),ZZ(0:LX),ZSQ(0:LX),Z2S2(0:LX)
       REAL (wp) :: ZS22(0:LX),X2(0:LX),ZZZ2(0:LX),ZSQ2(0:LX),ETA2(0:LX)
       REAL (wp) :: ZZ1(0:LX),ZZZ(0:LX)        
       REAL (wp) :: CX,ALC,TAX0,CLX,CRS,JTWALL,JZWALL 

!       COMMON /FUNX/XX,ITEILX                                           
!       COMMON /FUNX1/XX2,ITEILX1

       COMPLEX(wp) :: SUM1(0:LX),SUM2(0:LX),SUM3(0:LX),SUM4(0:LX)
       COMPLEX(wp) :: SUM5(0:LX),SUM6(0:LX),SUM7(0:LX),SUM8(0:LX)
       COMPLEX(wp) :: SUM9(0:LX),SUM10(0:LX),SUM11(0:LX),SUM12(0:LX)
       COMPLEX(wp) :: SUM13(0:LX),SUM14(0:LX),SUM15(0:LX),SUM16(0:LX)
       COMPLEX(wp) :: DSUM1(0:LX),DSUM2(0:LX),DSUM3(0:LX),SUM0(0:LX)

       REAL (wp) :: ETAA(0:LX),X(0:LX+2),ETAS(0:LX),RHO(0:LX)                    

       COMPLEX(wp) :: R1(0:LX),R2(0:LX),R3(0:LX)
       COMPLEX(wp) :: AAA(-1:LX),BBB(-1:LX)
       REAL (wp) :: XX(0:1000),XX2(0:1000)
       COMPLEX (wp) :: CFAK,CFAK1,S0,S1,S2,S3,S4,S5,S6

       INTEGER (wp) :: JANZ(-my:my,-160:160)
       REAL (wp) :: XT(ntp),YT(ntp),YZ(ntp), YF(ntp), YQ(ntp)
       COMPLEX (wp) :: I,VRALT,VZALT, VTIXJ,VZIXJ,VTIXJ1,VZIXJ1

       DATA I/(0.0d0,1.0d0)/
       LOGICAL IBAND,IBDEN
       REAL (wp) :: MU1(0:LX), MU2(0:LX)
       REAL (wp) :: NUE,XAA(15),MUEC,KAP(0:LX)

       COMPLEX (wp) :: AM1(2,2),AM2(2,2),AM3(2,2),AM4(2,2),AM5(2,2),            
     .         AV1(2),AV2(2),AV3(2)                                     
       COMPLEX (wp) :: AA(2,2,0:LX),BB(2,2,0:LX),CC(2,2,0:LX),
     .         DD(2,0:LX),EE(2,2,0:LX),FF(2,0:LX),UU(2,0:LX)                       
       COMPLEX (wp) :: ALJ(0:LX),BETJ(0:LX),GAMJ(0:LX),DELJ(0:LX)       
                                                                        
       REAL (wp) :: VEC(0:LX+1),XA1(10),XA2(10), KMN1, KMN2                     
                                                                        
       DATA XA1/0.0d0,12.0d0,12.0d0,0.0d0,0.0d0,2*6.0d0,
     .          18.0d0,18.0d0,6.0d0/
       DATA XA2/14.d0,26.0d0,26.0d0,14.d0,14.d0,2*6.0d0,
     .          18.0d0,18.0d0,6.0d0/ 
       DATA XAA/0.0d0,2*27.0d0,2*0.0d0,2*18.0d0,0.0d0,
     .          27.0d0,2*7.2d0,2*9.0d0,0.,18.d0/    
!                                                                       
!#! Marco, equilibrium parameters

       integer, parameter :: ipr=-2
       real (wp) , parameter :: chi=0.d0,beta0=0.d0
c theta=1.61 F=0.15
       real (wp) , parameter :: theta0 = 1.d0, alpha = 4.2d0
       real (wp) , parameter :: q01 = 0.d0, aa1 = 1.8748d0
       real (wp) , parameter :: bb1 = 0.83232d0
       real (wp) , parameter :: rr = 4.d0
                                                  
       integer :: ilinear,inorm,izeit
       DATA ILINEAR /0/, INORM /0/, IZEIT /0/                           
       REAL (wp) :: BZ0(-1:LX), BT0(-1:LX), P01(-1:LX)                                                                  
       INTEGER (wp) :: MPLO(3),NPLO(3),MPLOT(1)

       DATA MPLO/0,1,1/,NPLO/0,-3,-4/,MPLOT/3/                          

C         ITOR sottomultiplo di IB, altrimenti non funziona
c
       real(wp) :: xla, xle, c1
       real(wp) :: alfa1, alfa2, c11, c2, rs, r0
       real(wp) :: a, alf, gamma, epsi, bet, eps, eps1, beta, q0, b0
       real(wp) :: dt
       integer(wp) :: itend,ib,itor
       real (wp) :: eta0, alet, beet, rho0, gaet
       real (wp) :: mue, almu, bemu, kappa
       integer :: iacc, iband1, inopert, ibein, ibaus, ibout1

       DATA DT/0.01d0/
       DATA ITEND/10000/,IB/1000/,ITOR/2000000/

       DATA XLA/0.0d0/,XLE/1.0d0/,C1/0.d0/
                                                                        
       DATA ALFA1/4.0d0/,ALFA2/12.0d0/,C11/0.0d0/,C2/0.00d9/,
     .      RS/0.4d0/,R0/0.6d0/    
       DATA A/1.0d0/,ALF/0.6d0/,GAMMA/1.66666666d0/,EPSI/0.010d0/,
     .      BET/0.04d0/,
     .      EPS/0.33333333333333d0/,EPS1/1.0E-6/,
     .      BETA/1.d0/,Q0/0.7d0/,B0/1.0d0/
!                                                                       
! resistivita:                                                          
       DATA ETA0/3.3333E-5/,ALET/21.d0/,BEET/10.d0/,RHO0/1.0d0/,GAET/1.0d0/
! viscosita :                                                           
       DATA  MUE/3.3333E-2/,ALMU/1.d0/,BEMU/0.d0/, KAPPA/0.0E-6/
!                                                                       
       DATA IACC/ 1 /                                                   
!C      DATA IBAND/.TRUE./,IBAND1/0/,IBDEN/.TRUE./                       
!c inizia simulazione :
       DATA IBAND/.FALSE./,IBAND1/0/,INOPERT/0/
!c continua simulazione :
!       DATA IBAND/.TRUE./,IBAND1/0/,INOPERT/1/
                                                                        
!           LEGGE      1SOLO     TUTTO                                  
       DATA IBEIN/03/,IBAUS/04/,IBOUT1/02/                              
                                            
!#! Marco, other declarations
       integer(wp) :: iaus, ig
       real(wp) :: pi2, dt2, adt, dtmue,dtkap,gam1,p0
       integer(wp) :: itp,j,j0,n,im,m,nj,jmn,ms,ns0,ns1,j1
       integer(wp) :: zs,zdr,lxx,it,ns,ibib,ly,izz
       real(wp) :: one,dx,bzh,bth,p0wall
       real(wp) :: err,bz_avg,bzwall,btwall,a22,zeit
       real :: vrmax,vtmax,vzmax,brmax,btmax,bzmax
       real :: yvr,yvt,yvz,ybr,ybt,ybz,valmax
       real :: facnor, torf1,yfitp,yth,yqitp
       integer :: im0, im1,in0,in1,idx, IAM, js, k, l, ix
       real(wp) :: ampa, ampb,amp_old,phold,acontr,afluct
       real(wp) :: fak1,ampv,a12,alfs,bets,a11,a13,rrad
       real(wp) :: ftor, torflx, corr_fact, sign_r, sign_t
       real(wp) :: ybritp,yvritp,ybrfitp,yvrfitp,ybr1,ybr2,yvr1,yvr2
       real(wp) :: fsusi,ybtot,yvtot,ybf,yvf,torf
                                                                        
       integer :: i2d, inoconv, ivac, ippcd
       real(wp) :: lambda
       DATA I2D/0/, INOCONV/0/, IVAC/0/, IPPCD/0/
       DATA lambda/1.d0/
