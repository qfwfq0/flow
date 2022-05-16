       module defs

!-----------parameters-------------------------------------------

       integer , parameter :: wp = kind(1.d0) ! working precision
       integer (wp)  :: LX
       integer (wp)  :: my
       integer (wp) ,  parameter :: numz=600
       integer (wp) ,  parameter :: ntp=11000
       integer (wp)  :: nsec, nunitd, ios, ninput
       character(400) :: string, stringminus


!-----------scelta degli spettri---------------------------------

       integer (wp) :: n0 
       integer (wp) :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10
       integer (wp) :: n11, n12, n13, n14, n15, n16, n17, n18, n19, n20
       integer (wp) :: iz
       integer (wp) :: N2D 
       integer (wp) :: M2D 
       integer (wp) , allocatable :: MM(:),NANF(:),MWERT(:)
       integer (wp) , allocatable :: na1(:),na(:),ne(:)
       integer (wp) , allocatable :: nz(:),nzcon(:)
       INTEGER (wp) , allocatable :: CONVJ(:,:),CONVJ1(:,:)
       INTEGER (wp) , allocatable :: CONVMS(:,:),CONVNS(:,:)
       INTEGER (wp) , allocatable :: CONVNJ(:)
       INTEGER (wp) , allocatable :: JANZ(:,:)
       logical  :: fasi_sganciate
       
!-----------math-------------------------------------------

       real(wp) :: pi
       complex(wp) :: i
       integer(wp) :: seq ! to distinguish first or second application
!of Faraday law

!-----------I/O-------------------------------------------

       character(len=555) :: prefix, simdir
       character(len=555) :: spc2file, spc3file, ioiofile
       character(len=555) :: un_all, un_end, un_prevend
       character(len=555) :: unit_all, unit_end, unit_prevend
       character(len=555) :: spectrumfile
       CHARACTER(len=555) :: cwd
       integer(wp) :: dumint, dumintminus
       integer(wp) , parameter :: spectrumunit=5748


!-----------MP-------------------------------------------

       real(wp) :: dumr,dumt,dumz
       complex(wp) , allocatable :: bbr(:)
       complex(wp) , allocatable :: bbt(:), bbz(:)
       integer(wp) :: uft 
       integer(wp) :: statuss
       integer(wp) :: ik, loc, tempmv
       real(wp) :: temp_ph_mp
       logical :: file_exist
       character(400) :: ph_mp_file
!#! Marco, irrotational MP
       logical :: irrot
       real(wp) :: br_pert_freq
       integer(wp) :: qmp 
!#! amplitude and phase
       complex (wp) , allocatable :: mp_cmplx(:)
       real (wp) , allocatable  ::  mp_perc(:)
       real (wp) , allocatable  ::  ph_mp(:)
!#! m_MP and n_MP (2D case)
       integer(wp) , allocatable :: m_mp(:)
       integer(wp) , allocatable :: n_mp(:)
       integer(wp) , allocatable :: j_mp(:), temp_mp(:)
       real (wp) , allocatable :: phmp(:), amp_a(:)
!#! old MP
       integer(wp) , allocatable :: m_mpold(:)
       integer(wp) , allocatable :: n_mpold(:) 
       integer(wp) , allocatable :: j_mpold(:)
       real (wp) , allocatable :: ampold(:)
       real (wp) :: pitol
!#!
       integer(wp) :: mloc, nloc 
!#! Marco, rotational MP for tokamak simulations
       logical :: tok_rot,mp_exp,mp_tanh
       integer(wp) :: tkp
       real(wp) :: tau_growth,tau_decrease,time_passage
       real(wp) :: time_medium
       integer(wp) :: tot_it, it_beginning
    
!-----------Momentum source-------------------------------------------

       real(wp) :: m0_momsour(0:2)
       integer(wp) :: momsour_shape
       real(wp) :: gauss_maxvz0,gauss_fwhm,gauss_center
       complex(wp) , allocatable :: momsour(:,:)

!-----------Dump velocity--------------------------------------------
       real(wp) :: dump_v 

!-----------arrays declaration-------------------------------------------

!#! Marco, for backward compatibility
       real , allocatable :: BROLD(:,:),BTOLD(:,:)
     .             ,BZOLD(:,:)   
       real , allocatable :: VROLD(:,:),VTOLD(:,:)
     .             ,VZOLD(:,:)   
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
       COMPLEX (wp) , allocatable :: AT(:,:),AZ(:,:)
       COMPLEX (wp) , allocatable :: ATOLD(:,:),AZOLD(:,:)
       
       COMPLEX (wp) :: bzaini,bzafin
       REAL (wp) :: E0,VRLX,B2
       REAL (wp) , allocatable :: BREQ(:), BTEQ(:), BZEQ(:)   
       REAL (wp) , allocatable :: VREQ(:), VTEQ(:), VZEQ(:)   
       REAL (wp) , allocatable :: ohmbt(:), ohmbz(:)
       REAL (wp) , allocatable :: ohmbt_old(:), ohmbz_old(:)
       REAL (wp) , allocatable :: dummy(:)
       REAL (wp) , allocatable :: evr0(:),ebt0(:),ebz0(:)
       REAL (wp) , allocatable :: vr0(:),vz0(:)
       REAL (wp) , allocatable :: ZS2(:),ZZ(:),ZSQ(:),Z2S2(:)
       REAL (wp) , allocatable :: ZS22(:),X2(:),ZZZ2(:),ZSQ2(:),ETA2(:)
       REAL (wp) , allocatable :: ZZ1(:),ZZZ(:)        
       
       REAL (wp) :: CX,ALC,TAX0,CLX,CRS,JTWALL,JZWALL 

       COMPLEX(wp) , allocatable :: SUM1(:),SUM2(:),SUM3(:),SUM4(:)
       COMPLEX(wp) , allocatable :: SUM5(:),SUM6(:),SUM7(:),SUM8(:)
       COMPLEX(wp) , allocatable :: SUM9(:),SUM10(:),SUM11(:),SUM12(:)
       COMPLEX(wp) , allocatable :: SUM13(:),SUM14(:),SUM15(:),SUM16(:)
       COMPLEX(wp) , allocatable :: DSUM1(:),DSUM2(:),DSUM3(:),SUM0(:)

       REAL (wp) , allocatable :: ETAA(:),X(:),ETAS(:),RHO(:) 
       REAL (wp) , allocatable :: XX(:),XX2(:)
       COMPLEX(wp) , allocatable :: R1(:),R2(:),R3(:)
       COMPLEX(wp) , allocatable :: AAA(:),BBB(:)
       REAL (wp) , allocatable :: MU1(:), MU2(:),KAP(:)
       COMPLEX (wp) :: CFAK,CFAK1,S0,S1,S2,S3,S4,S5,S6

       REAL (wp) :: XT(ntp),YT(ntp),YZ(ntp), YF(ntp), YQ(ntp)
       COMPLEX (wp) :: VRALT,VZALT, VTIXJ,VZIXJ,VTIXJ1,VZIXJ1

       LOGICAL :: IBAND
       REAL (wp) :: NUE,XAA(15),MUEC
       COMPLEX (wp) :: AM1(2,2),AM2(2,2),AM3(2,2),AM4(2,2),AM5(2,2),            
     .         AV1(2),AV2(2),AV3(2)                                     
       COMPLEX (wp) , allocatable :: AA(:,:,:),BB(:,:,:),CC(:,:,:),
     .         DD(:,:),EE(:,:,:),FF(:,:),UU(:,:)                       
       COMPLEX (wp) , allocatable :: ALJ(:),BETJ(:),GAMJ(:),DELJ(:)       
                                                                        
       REAL (wp) , allocatable :: VEC(:)
       REAL (wp) :: XA1(10),XA2(10), KMN1, KMN2                     

!-----------equilibrium-------------------------------------------

       integer(wp) :: ipr
       real (wp)  :: chi, beta0
       real (wp)  :: theta0, alpha, e0_over_eta, ohmbz0
       real (wp)  :: q01, aa1
       real (wp)  :: bb1
       real (wp)  :: rr
       integer :: ilinear,inorm,izeit
       REAL (wp) , allocatable :: BZ0(:), BT0(:), P01(:)
       INTEGER (wp) :: MPLO(3),NPLO(3),MPLOT(1)
       real(wp) :: xla, xle, c1
       real(wp) :: alc1, alc2, alfa1, alfa2, c11, c2, rs, r0
       real(wp) :: a, alf, gamma, epsi, bet, eps, eps1, beta, q0, b0
       real(wp) :: dt, cx1, cx2
       integer(wp) :: itend,ib,itor
       real (wp) :: eta0, alet, beet, rho0, gaet
       real (wp) :: mue, almu, bemu, kappa
       integer :: iacc, iband1, inopert

!#! Marco, importantissimi per output, non cambiare :)
       integer(wp) , parameter :: ibein = 03
       integer(wp) , parameter :: ibaus = 04 
       integer(wp) , parameter :: ibout1 = 02

!#! Marco, other declarations
       integer(wp) :: iaus, ig
       integer(wp) :: nx1, ixxx
       real(wp) :: pi2, dt2, adt, dtmue,dtkap,gam1,p0, xst
       integer(wp) :: itp,j,j0,n,im,m,nj,jmn,ms,ns0,ns1,j1
       integer(wp) :: zs,zdr,lxx,it,ns,ibib,z2s
       integer(wp) :: ly, izz
       real(wp) :: one,dx,bzh,bth,p0wall,xxxx,xxx1,dr
       REAL :: HELI2D
       real(wp) :: err,bz_avg,bzwall,btwall,a22,zeit
!#! Marco, non li dichiaro doubles perché sennò casino
       real :: vrmax,vtmax,vzmax,brmax,btmax,bzmax
       real :: yvr,yvt,yvz,ybr,ybt,ybz,valmax
       real :: facnor, torf1,yfitp,yth,yqitp
       integer :: im0, im1,in0,in1,idx, IAM, js, k, l, ix
       real(wp) :: ampa, ampb,amp_old,phold,acontr,afluct
       real(wp) :: fak1,ampv,a12,alfs,bets,a11,a13,rrad
       real(wp) :: ftor, torflx, corr_fact, sign_r, sign_t
       real(wp) :: ybritp,yvritp,ybrfitp,yvrfitp,ybr1,ybr2,yvr1,yvr2
       real(wp) :: fsusi,ybtot,yvtot,ybf,yvf,torf
       integer(wp) :: i2d, inoconv, ivac, ippcd
       real(wp) :: lambda

       contains 
       subroutine alloc_arrays(amy,aiz,alx)
   
        implicit none

        integer(wp), intent(in) :: amy, aiz, alx

        allocate( MM(0:amy), NANF(0:amy), MWERT(-amy:amy) )
        allocate( na1(0:amy), na(-amy:amy), ne(-amy:amy) )
        allocate( nz(0:amy), nzcon(0:amy) )
        allocate( CONVJ(0:aiz,0:2*aiz+1), CONVJ1(0:aiz,0:2*aiz+1) )
        allocate( CONVMS(0:aiz,0:2*aiz+1), CONVNS(0:aiz,0:2*aiz+1) )
        allocate( CONVNJ(0:aiz) )
        allocate( JANZ(-amy:amy,-900:900) )
        allocate( bbr(0:alx) )
        allocate( bbt(-1:alx), bbz(-1:alx) )
        allocate( mp_perc(0:qmp), mp_cmplx(0:qmp) )
        allocate( ph_mp(0:qmp) )
        allocate( m_mp(0:qmp) )
        allocate( n_mp(0:qmp) )
        allocate( j_mp(0:qmp), temp_mp(0:qmp) )
        allocate( phmp(0:qmp), amp_a(0:qmp) )
        allocate( m_mpold(0:qmp) )
        allocate( n_mpold(0:qmp) )
        allocate( j_mpold(0:qmp) )
        allocate( ampold(0:qmp) )
        allocate( momsour(0:alx,0:2) )
!#! Marco, for backward compatibility
        allocate( BROLD(0:alx,-aiz:aiz),BTOLD(-1:alx,-aiz:aiz)
     .           ,BZOLD(-1:alx,-aiz:aiz))
        allocate( VROLD(0:alx,-aiz:aiz),VTOLD(-1:alx,-aiz:aiz)
     .           ,VZOLD(-1:alx,-aiz:aiz))  
        allocate( BR(0:alx,-aiz:aiz),BT(-1:alx,-aiz:aiz)
     .           ,BZ(-1:alx,-aiz:aiz))
        allocate( VR(0:alx,-aiz:aiz),VT(-1:alx,-aiz:aiz)
     .           ,VZ(-1:alx,-aiz:aiz))  
        allocate( BRN(0:alx,-aiz:aiz),BTN(-1:alx,-aiz:aiz)
     .           ,BZN(-1:alx,-aiz:aiz))   
        allocate( VRN(0:alx,-aiz:aiz),VTN(-1:alx,-aiz:aiz)
     .           ,VZN(-1:alx,-aiz:aiz))
        allocate( VRP(0:alx,-aiz:aiz),VTP(-1:alx,-aiz:aiz))
        allocate( JR(-1:alx,-aiz:aiz),JT(0:alx,-aiz:aiz)
     .           ,JZ(0:alx,-aiz:aiz))
        allocate( ER(0:alx,-aiz:aiz), ET(0:alx,-aiz:aiz)
     .           ,EZ(0:alx,-aiz:aiz))
        allocate( PHI(0:alx,-aiz:aiz),RHOC(0:alx,-aiz:aiz))
        allocate( AT(-1:alx,-aiz:aiz),AZ(-1:alx,-aiz:aiz))
        allocate(ATOLD(-1:alx,-aiz:aiz),AZOLD(-1:alx,-aiz:aiz))
        allocate( BREQ(0:alx), BTEQ(-1:alx), BZEQ(-1:alx) )
        allocate( VREQ(0:alx), VTEQ(-1:alx), VZEQ(-1:alx) )   
        allocate( ohmbt(0:alx+1), ohmbz(0:alx+1) )
        allocate( ohmbt_old(0:alx+1), ohmbz_old(0:alx+1) )
        allocate( dummy(0:alx+1) )
        allocate( evr0(0:alx), ebt0(-1:alx), ebz0(-1:alx) )
        allocate( vr0(0:alx), vz0(0:alx) )
        allocate( ZS2(0:alx), ZZ(0:alx), ZSQ(0:alx), Z2S2(0:alx) )
        allocate( ZS22(0:alx), X2(0:alx), ZZZ2(0:alx) )
        allocate( ZSQ2(0:alx), ETA2(0:alx) )
        allocate( ZZ1(0:alx),ZZZ(0:alx) )        
        allocate( SUM1(0:alx), SUM2(0:alx), SUM3(0:alx), SUM4(0:alx) )
        allocate( SUM5(0:alx), SUM6(0:alx), SUM7(0:alx), SUM8(0:alx) )
        allocate( SUM9(0:alx), SUM10(0:alx), SUM11(0:alx) )
        allocate( SUM13(0:alx), SUM14(0:alx), SUM15(0:alx) )
        allocate( DSUM1(0:alx), DSUM2(0:alx), DSUM3(0:alx))
        allocate( SUM12(0:alx) , SUM16(0:alx), SUM0(0:alx) )
        allocate( ETAA(0:alx), X(0:alx+2), ETAS(0:alx), RHO(0:alx) ) 
        allocate( XX(0:1000), XX2(0:1000) )
        allocate( R1(0:alx), R2(0:alx), R3(0:alx) )
        allocate( AAA(-1:alx), BBB(-1:alx) )
        allocate( MU1(0:alx), MU2(0:alx), KAP(0:alx) )
        allocate( AA(2,2,0:alx), BB(2,2,0:alx), CC(2,2,0:alx) )
        allocate( DD(2,0:alx), EE(2,2,0:alx), FF(2,0:alx), UU(2,0:alx) )
        allocate( ALJ(0:alx), BETJ(0:alx), GAMJ(0:alx), DELJ(0:alx) )       
        allocate( VEC(0:alx+1) )
        allocate( BZ0(-1:alx), BT0(-1:alx), P01(-1:alx) )

       end subroutine alloc_arrays

       subroutine dealloc_arrays
   
        implicit none

        deallocate( MM, NANF, MWERT )
        deallocate( na1, na, ne )
        deallocate( nz, nzcon )
        deallocate( CONVJ, CONVJ1 )
        deallocate( CONVMS, CONVNS )
        deallocate( CONVNJ )
        deallocate( JANZ )
        deallocate( bbr )
        deallocate( bbt, bbz )
        deallocate( mp_perc, mp_cmplx )
        deallocate( ph_mp )
        deallocate( m_mp )
        deallocate( n_mp )
        deallocate( j_mp, temp_mp )
        deallocate( phmp, amp_a )
        deallocate( m_mpold )
        deallocate( n_mpold )
        deallocate( j_mpold )
        deallocate( ampold )
        deallocate( momsour )
        deallocate( BR, BT, BZ )
        deallocate( VR, VT, VZ )  
!#! Marco, for backward compatibility
        deallocate( BROLD, BTOLD, BZOLD )
        deallocate( VROLD, VTOLD, VZOLD )  
        deallocate( BRN, BTN, BZN )   
        deallocate( VRN, VTN, VZN )
        deallocate( VRP, VTP )
        deallocate( JR,JT, JZ )
        deallocate( ER, ET, EZ )
        deallocate( PHI, RHOC )
        deallocate( AT, AZ )
        deallocate( ATOLD, AZOLD )
        deallocate( BREQ, BTEQ, BZEQ )
        deallocate( VREQ, VTEQ, VZEQ )   
        deallocate( ohmbt, ohmbz )
        deallocate( ohmbt_old, ohmbz_old )
        deallocate( dummy )
        deallocate( evr0, ebt0, ebz0 )
        deallocate( vr0, vz0 )
        deallocate( ZS2, ZZ, ZSQ, Z2S2 )
        deallocate( ZS22, X2, ZZZ2, ZSQ2, ETA2 )
        deallocate( ZZ1,ZZZ )        
        deallocate( SUM1, SUM2, SUM3, SUM4 )
        deallocate( SUM5, SUM6, SUM7, SUM8 )
        deallocate( SUM9, SUM10, SUM11, SUM12 )
        deallocate( SUM13, SUM14, SUM15, SUM16 )
        deallocate( DSUM1, DSUM2, DSUM3, SUM0 )
        deallocate( ETAA, X, ETAS, RHO ) 
        deallocate( XX, XX2 )
        deallocate( R1, R2, R3 )
        deallocate( AAA, BBB )
        deallocate( MU1, MU2, KAP )
        deallocate( AA, BB, CC )
        deallocate( DD, EE, FF, UU )
        deallocate( ALJ, BETJ, GAMJ, DELJ )       
        deallocate( VEC )
        deallocate( BZ0, BT0, P01 )

       end subroutine dealloc_arrays
  
       subroutine init_variab

!#! Marco, initialize necessary variables (some of them were never
!changed in the past, and I don't understand their meaning)

         pi = 3.14159265358979323846d0
         i = (0.d0,1.d0)
         seq = 1
         XA1 = (/ 0.0d0,12.0d0,12.0d0,0.0d0,0.0d0,6.0d0,
     .          6.0d0,18.0d0,18.0d0,6.0d0 /)
         XA2 = (/ 14.d0,26.0d0,26.0d0,14.d0,14.d0,6.0d0,
     .          6.0d0,18.0d0,18.0d0,6.0d0 /)
         XAA = (/ 0.0d0,27.0d0,27.0d0,0.d0,0.0d0,18.0d0,
     .           18.0d0,0.0d0,27.0d0,7.2d0,7.2d0,9.0d0,9.0d0,0.d0
     .          ,18.d0 /)

         ILINEAR = 0 ; INORM = 0 ; IZEIT = 0 
         MPLO = (/ 0,1,1 /) ; NPLO = (/ 0,-3,-4 /) 
         MPLOT = 3

         XLA = 0.0d0 ; XLE = 1.0d0 ; C1 = 0.d0
         ALFA1 = 4.0d0 ; ALFA2 = 12.0d0 ; C11 = 0.0d0  
         C2 = 0.0d0 ; RS = 0.4d0 ; R0 = 0.6d0
         A = 1.0d0 ; ALF = 0.6d0 ; GAMMA = 1.66666666d0
         EPSI = 0.010d0 ; BET = 0.04d0 ; EPS = 0.33333333333333d0
         EPS1 = 1.0E-6 ; BETA = 1.d0 ; Q0 = 0.7d0 ; B0 = 1.0d0

         iacc = 1 
         INOCONV = 0 ; IVAC = 0 ;  IPPCD = 0
         lambda = 1.d0

         do k=0,qmp
!#! Marco, inizializzo così, tanto poi viene ricalcolato
          j_mp(k) = -1
          j_mpold(k) = -1
         enddo
         mloc = 0
         nloc = 0
       end subroutine init_variab
 
       subroutine init_arrays(amy,aiz,alx)
   
        implicit none

        integer(wp), intent(in) :: amy, aiz, alx

!#! Marco, initialize big arrays

        br(0:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        bt(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        bz(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        vr(0:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        vt(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        vz(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        brn(0:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        btn(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        bzn(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        vrn(0:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        vtn(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
        vzn(-1:alx,-aiz:aiz) = dcmplx(0.d0, 0.d0) 
!#! Marco, for backward compatibility
        brold(0:alx,-aiz:aiz) = 0.d0 
        btold(-1:alx,-aiz:aiz) =0.d0 
        bzold(-1:alx,-aiz:aiz) =0.d0 
        vrold(0:alx,-aiz:aiz) = 0.d0 
        vtold(-1:alx,-aiz:aiz) =0.d0 
        vzold(-1:alx,-aiz:aiz) =0.d0 
        er(:,:) = dcmplx(0.d0, 0.d0) 
        et(:,:) = dcmplx(0.d0, 0.d0) 
        ez(:,:) = dcmplx(0.d0, 0.d0) 
        jr(:,:) = dcmplx(0.d0, 0.d0) 
        jt(:,:) = dcmplx(0.d0, 0.d0) 
        jz(:,:) = dcmplx(0.d0, 0.d0) 
        vrp(:,:) = dcmplx(0.d0, 0.d0) 
        vtp(:,:) = dcmplx(0.d0, 0.d0) 
        phi(:,:) = dcmplx(0.d0, 0.d0) 
        rhoc(:,:) = dcmplx(0.d0, 0.d0) 
        at(:,:) = dcmplx(0.d0, 0.d0) 
        az(:,:) = dcmplx(0.d0, 0.d0) 
        atold(:,:) = dcmplx(0.d0, 0.d0) 
        azold(:,:) = dcmplx(0.d0, 0.d0) 

       end subroutine init_arrays

       subroutine default_spc2 
  
!#! Marco, default parameters, in case spc2 is not read well
        lx = 100
        my = 10 
        iz = 10 
        i2d = 1 
        n2d = -10 
        m2d = 1 
        fasi_sganciate = .false.
        irrot = .true.
       end subroutine default_spc2 

       subroutine default_ioio
        
        prefix="/0/rfx/ricercatori/ft/specyl/veranda/flow" 
        simdir="/archive/rfp/3d/theta16/base/" 
        un_all="specyl_1_all.dat"
        un_end="specyl_1_end.dat"
        un_prevend='dum_0_end.dat'

       end subroutine default_ioio

       subroutine default_spc3(mi2d)
  
        integer(wp), intent(in) :: mi2d
!#! Marco, default parameters, in case spc3 is not read well
!        if (mi2d .eq. 1) then 
!!        mm(1) = 1
!!        mm(2) = 2
!!        mm(3) = 3
!!        mm(4) = 4
!!        mm(5) = 5
!!        mm(6) = 6
!!        mm(7) = 7
!!        mm(8) = 8
!!        mm(9) = 9
!!        mm(10) = 10
!         mm = (/ 0,1,2,3,4,5,6,7,8,9,10 /)
!         nanf = (/ 0,-10,-20,-30,-40,-50,-60,-70,-80,-90,-100 /)
!         nz = (/ 0,1,1,1,1,1,1,1,1,1,1 /)
!         nzcon = (/ 1,1,1,1,1,1,1,1,1,1,1 /)
!        else
!!         do k=0,mmy 
!!          mm(k) = k
!!         enddo 
!         mm = (/ 0,1,2,3,4,5,6,7,8,9,10 /)
!         nanf = (/ -25,-55,-50,-60,-70,-80,-90,-100,-110,-120,-130 /)
!         nz = (/ 25,65,45,45,45,45,45,45,45,45,45 /)
!         nzcon = (/ 26,65,45,45,45,45,45,45,45,45,45 /)
!        endif

        n0 = 0
        n1 = 1
        n2 = 0
        n3 = 0
        n4 = 0
        n5 = 0
        n6 = 0
        n7 = 0
        n8 = 0
        n9 = 0
        n10 = 0
        n11 = 0
        n12 = 0
        n13 = 0
        n14 = 0
        n15 = 0
        n16 = 0
        n17 = 0
        n18 = 0
        n19 = 0
        n20 = 0
        qmp = 0
        br_pert_freq = 0.d0
        mp_cmplx = cmplx(0.d0,0.d0)
        mp_perc = 0.d0
        ph_mp = 0.d0
        m_mp = 0
        n_mp = 0
        m_mpold = 0 
        n_mpold = 0 
        ampold = 0.d0
        pitol = 1.e-4
        tok_rot = .false.
        mp_exp = .true.
        mp_tanh = .false.
        tkp = 29
        tot_it = 0
        tau_growth = 2.d-10 !very slow growth time scale means no effect
! on a fixed MP
        tau_decrease = 2.d-10
        time_passage = 1.d10 !very high passage time means the field
        time_medium = 2.d10 !very high passage time means the field
!will never start decreasing
        m0_momsour(0) = 0.d0
        m0_momsour(1) = 0.d0
        m0_momsour(2) = 0.d0
!#! momsour_shape:
!#! momsour_shape = 0 : constant, in z direction only for now
!#! momsour_shape = 1 : it gives a gaussian vz profile
        momsour = 0
!#! in case of momsour_shape=1
        gauss_maxvz0 = 0.02
        gauss_fwhm = 0.05
        gauss_center = 0.50
!#!dump velocity harmonic; if dump_v=1. do nothing. If dump_v. ne. 1
!dump velocity of the factor (0<dump_v<infty)
        dump_v = 1.
!#! other stuff
        ipr = -2
        chi = 0.d0
        beta0 = 0.d0
        theta0 = 1.d0
        alpha = 4.d0
        q01 = 0.d0
        aa1 = 1.8748d0
        bb1 = 0.83232d0
        rr = 4.d0
        dt = 0.01d0
        itend = 10000
        ib = 100
        itor = 1000000
        eta0 = 3.333e-5
        alet = 21.d0
        beet = 10.d0
        rho0 = 1.d0
        gaet = 1.d0
        mue = 3.333e-2
        almu = 1.d0
        bemu = 0.d0
        lambda = 0.d0
        iband = .false. ! inizia la simulazione
!        iband = .true. ! continua la simulazione
        iband1 = 0
        inopert = 0 ! perturba la simulazione
!        inopert = 1 ! non perturba la simulazione

       end subroutine default_spc3

c     int2char
c     ################################################################
      function int2char(n) result (chr)

      implicit none

      integer(8)      :: n
      character(10):: chr

      integer      :: i,exp,k,j
      character(3) :: c

      if (abs(n) > 0) then
         exp = int(log(float(n))/log(1d1))
      else
         exp = 0
      endif

      if (n >= 0) then
        chr=''
      else
        chr='-'
      endif

      k = abs(n)
      do i=exp,0,-1
         j = k/10**i
         c = achar(48+j)
         chr = trim(chr)//trim(c)
         if (i > 0 .and. j /= 0) k = mod(k,j*10**i)
      enddo

      end function int2char

      subroutine str2int(str,int,stat)
          implicit none
          ! Arguments
          character(len=*),intent(in) :: str
          integer(wp),intent(out)         :: int
          integer(wp),intent(out)         :: stat
      
          read(str,*,iostat=stat)  int
       end subroutine str2int
      
       end module defs
