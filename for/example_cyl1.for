C  MY=numero degli M(senza NULL)                
C  N0=numero degli N (M=0) (senza NULL)                                 
C  N1...N30= numero degli  N con tutti gli  MY (con il NULL)            
C  NANF = N-ANFANG, BEI NANF(0) OHNE NULL, SONST MIT NULL               
C   NZ(30), GROESSTES M AUF 100 BEGRENZT                                
                                                                        
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
       DATA   PI/3.14159265358979323846d0/
 
       real(wp) t(ntp)                 
                                               
       INTEGER  :: MM(0:MY),NANF(0:MY),MWERT(-my:my)
       INTEGER  :: NA1(0:my),NA(-my:my),NE(-my:my)                                                  
!       DATA MM/0,1,2,3,4/
!       DATA NANF/-25,-55,-50,-60,-70/ 
       DATA MM/0,1,2/
       DATA NANF/-25,-55,-50/ 

!#! Marco, perturbazione iniziale con fasi sganciate
       logical ,parameter :: fasi_sganciate=.false.
!#! Marco, irrotational MP
       logical ,parameter :: irrot=.true.
       real(wp), parameter :: br_pert_freq=0.d0
       integer, parameter :: qmp = 0
!#! amplitude and phase
       real (wp) , dimension(0:qmp), parameter ::
     .             mp_perc=(/6.d-3/)
       real (wp) , dimension(0:qmp), parameter :: 
     .             ph_mp=(/0.d0/)
!#! m_MP and n_MP (2D case)
       integer , parameter :: m_mp = 1
       integer , dimension(0:qmp), parameter :: n_mp = (/-8/)
       integer , dimension(0:qmp) :: j_mp, temp_mp
       real (wp) , dimension(0:qmp) :: phmp, amp_a
!#! old MP
       integer , parameter :: m_mpold = 0
       integer , dimension(0:qmp), parameter :: n_mpold = (/0/)
       integer , dimension(0:qmp) :: j_mpold
       real (wp) , dimension(0:qmp), parameter :: 
     .             ampold=(/0.d0/)
       real (wp) :: pitol
       data   pitol/1.e-5/
    
!#! Marco, momentum source
       real(wp) :: m0_momsour(0:2)
       data m0_momsour/0.d0,0.d0,0.d0/
       complex(wp) ::momsour(0:lx,0:2)

!#! Marco, trick to break phase constanc
       logical , parameter :: jsa=.false.
       real (wp) ::  corr
       integer (wp) :: jstar
       COMPLEX(wp) :: BR(0:LX,-IZ:IZ),BT(-1:LX,-IZ:IZ)
     .             ,BZ(-1:LX,-IZ:IZ)   
       COMPLEX(wp) :: VR(0:LX,-IZ:IZ),VT(-1:LX,-IZ:IZ)
     .             ,VZ(-1:LX,-IZ:IZ)   
       COMPLEX(wp) :: BRN(0:LX,-IZ:IZ),BTN(-1:LX,-IZ:IZ)
     .             ,BZN(-1:LX,-IZ:IZ)   
       COMPLEX(wp) :: VRN(0:LX,-IZ:IZ),VTN(-1:LX,-IZ:IZ)
     .             ,VZN(-1:LX,-IZ:IZ) 
       COMPLEX(wp) :: VRP(0:LX,-IZ:IZ),VTP(-1:LX,-IZ:IZ)                      
       COMPLEX(wp) :: JR(-1:LX,-IZ:IZ),JT(0:LX,-IZ:IZ)
     .             ,JZ(0:LX,-IZ:IZ)   
       COMPLEX(wp) :: ER(0:LX,-IZ:IZ), ET(0:LX,-IZ:IZ)
     .             ,EZ(0:LX,-IZ:IZ)
       COMPLEX(wp) :: PHI(0:LX,-IZ:IZ),RHOC(0:LX,-IZ:IZ)  
! 27/11/03 Daniele: per il calcolo di Az e Atheta (nel gauge Ar=0)
	COMPLEX (wp) :: AT(-1:LX,-IZ:IZ),AZ(-1:LX,-IZ:IZ)
! 27/02/04 Daniele: per il calcolo di Az e Atheta vecchi
	COMPLEX (wp) :: ATOLD(-1:LX,-IZ:IZ),AZOLD(-1:LX,-IZ:IZ)

       REAL (wp) :: ZS2(0:LX),ZZ(0:LX),ZSQ(0:LX),Z2S2(0:LX)                     
       REAL (wp) :: ZS22(0:LX),X2(0:LX),ZZZ2(0:LX),ZSQ2(0:LX),ETA2(0:LX)
       REAL (wp) :: ZZ1(0:LX),ZZZ(0:LX)                                         
       REAL (wp) :: CX,ALC,TAX0,CLX,CRS,JTWALL,JZWALL                           
                                                                        
       COMMON /FCTX/CX1,CX2,ALC1,ALC2,TAX0,IXXX,NX1,CLX,CRS             
       COMMON /FUNX/XX,ITEILX                                           
       COMMON /FUNX1/XX2,ITEILX1

       COMPLEX (wp) :: SUM1(0:LX),SUM2(0:LX),SUM3(0:LX),SUM4(0:LX),
     . SUM5(0:LX),
     . SUM6(0:LX),SUM7(0:LX),SUM8(0:LX),SUM9(0:LX),SUM10(0:LX),         
     . SUM11(0:LX),SUM12(0:LX),SUM13(0:LX),SUM14(0:LX),                 
     . SUM15(0:LX),SUM16(0:LX),                                         
     . DSUM1(0:LX),DSUM2(0:LX),DSUM3(0:LX),SUM0(0:LX)                   
                                                                        
       REAL (wp) :: ETAA(0:LX),X(0:LX+2),ETAS(0:LX),RHO(0:LX)                    
                                                                        
       COMPLEX (wp) :: R1(0:LX),R2(0:LX),R3(0:LX),AAA(-1:LX),BBB(-1:LX) 
                                                                        
       REAL (wp) :: XX(0:1000),XX2(0:1000)                                      
                                                                        
       COMPLEX (wp) :: CFAK,CFAK1,S0,S1,S2,S3,S4,S5,S6                          
                                                                        
       CHARACTER TEXT*17,TEXT1*16,TEXT2*11,TEXT3*11,TEXTA*25            
                                                                        
       INTEGER (wp) :: JANZ(-my:my,-160:160),NZ(0:my),NZCON(0:my)

       REAL (wp) :: XT(ntp),YT(ntp),YZ(ntp), YF(ntp), YQ(ntp)                   
                                                                        
       COMPLEX (wp) :: I,VRALT,VZALT, VTIXJ,VZIXJ,VTIXJ1,VZIXJ1                 
                                                                        
       DATA I/(0.0d0,1.0d0)/                                                
                                                                        
       LOGICAL IBAND,IBDEN                                              
                                                                        
       REAL (wp) :: MUE, ALMU, BEMU, MU1(0:LX), MU2(0:LX)                       
       REAL (wp) :: NUE,XAA(15),MUEC,KAPPA,KAP(0:LX)                  

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
                                                  

       DATA ILINEAR /0/, INORM /0/, IZEIT /0/                           
       REAL (wp) :: BZ0(-1:LX), BT0(-1:LX), P01(-1:LX)                          
!                                                                       
       INTEGER (wp) :: MPLO(3),NPLO(3)                                          
                                                                        
       DATA MPLO/0,1,1/,NPLO/0,-3,-4/,MPLOT/3/                          
                                                                        
C         ITOR sottomultiplo di IB, altrimenti non funziona
c
       DATA DT/0.01/
       DATA ITEND/500000/,IB/1000/,ITOR/2000000/

       DATA XLA/0.0d0/,XLE/1.0d0/,C1/0.d0/
                                                                        
       DATA ALFA1/4.0d0/,ALFA2/12.0d0/,C11/0.0d0/,C2/0.00d9/,
     .      RS/0.4d0/,R0/0.6d0/    
       DATA A/1.0d0/,ALF/0.6d0/,GAMMA/1.66666666d0/,EPSI/0.010d0/,
     .      BET/0.04d0/,
     .      EPS/0.33333333333333d0/,EPS1/1.0E-6/,
     .      BETA/1.d0/,Q0/0.7d0/,B0/1.0d0/
!                                                                       
! resistivita:                                                          
       DATA ETA0/1.0E-5/,ALET/21.d0/,BEET/10.d0/,RHO0/1.0d0/,GAET/1.0d0/
! viscosita :                                                           
       DATA  MUE/1.0E-2/,ALMU/1.d0/,BEMU/0.d0/, KAPPA/0.0E-6/
!                                                                       
       DATA IACC/ 1 /                                                   
C      DATA IBAND/.TRUE./,IBAND1/0/,IBDEN/.TRUE./                       
c inizia simulazione :
       DATA IBAND/.FALSE./,IBAND1/0/,INOPERT/0/
c continua simulazione :
!       DATA IBAND/.TRUE./,IBAND1/0/,INOPERT/1/
                                                                        
!           LEGGE      1SOLO     TUTTO                                  
       DATA IBEIN/03/,IBAUS/04/,IBOUT1/02/                              
                                            
                                                                        
       DATA TEXTA/'( CYL *MHD )$'/

       DATA I2D/0/
       DATA INOCONV/0/
       DATA IVAC/0/
       DATA IPPCD/0/
       DATA lambda/1./
