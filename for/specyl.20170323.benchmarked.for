      PROGRAM SPECYL

      use defs
      use mp
      use mom_sour
      use math
      use equilibrium

      implicit none
      external getarg
      integer :: omp_get_num_threads,omp_get_thread_num

      if ( iargc() .eq. 1 ) then
         call getarg(1, string); read(string,*) nsec
         write(*,*) 'section number : #',trim(string)
      else
         write(*,*) "usage is: "
      end if

      namelist /spc2/ lx,my,iz,i2d,n2d,m2d
     $   ,fasi_sganciate,irrot

!#! I define a new namelist containing relevant infos about spectra,
!equilibrium, etc
      namelist /spc3/ qmp,br_pert_freq,mp_perc,ph_mp,m_mp,n_mp
     $   ,mm,nanf,nz,nzcon,m0_momsour
     $   ,m_mpold,n_mpold,ampold,pitol,ipr
     $   ,chi,beta0,theta0,alpha,q01,aa1,bb1,rr,dt
     $   ,itend,ib,itor,eta0,alet,beet,rho0,gaet,mue,almu,bemu
     $   ,kappa,iband,iband1,inopert
     $   ,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10

!#! preliminary I-O operations
      call getcwd(prefix)
!      write(*,*) 'prefix=',trim(prefix)
      unit_all = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_all.dat")
      unit_end = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_end.dat")
!      write(*,*) 'debug_files_0 ',trim(unit_all)
!      write(*,*) 'debug_files_0 ',trim(unit_end)

!#! Marco, load default values of the spc2 namelist
      call default_spc2 
 
      spc2file = trim(trim(prefix)//"/for/"//
     .      trim(string)//"/spc2.in")
      spc3file = trim(trim(prefix)//"/for/"//
     .      trim(string)//"/spc3.in")

!      write(*,*) 'spc2in',spc2file
      open(unit=nunitd, file=trim(spc2file)
     .     ,status='unknown', iostat=ios)
      if (ios /= 0) then
       write(*,*) 'could not open namelist spc2.in'
       stop
      endif
      read(nunitd,spc2,iostat=ios)
      if (ios /= 0) then
       write(*,*) 'could not read namelist spc2.in'
      else
       close(nunitd)
      endif

!#! Marco, after reading the fundamental parameters (my,lx,iz) 
!I do allocate the arrays
      call alloc_arrays(my,iz,lx)
!#! and then load the default values of many variables
      call default_spc3(i2d)

!#! I read the spc3 namelist
      nunitd=22
      open(unit=nunitd, file=trim(spc3file)
     .     ,status='unknown', iostat=ios)
      if (ios /= 0) then
       write(*,*) 'could not open namelist spc3.in'
       stop
      endif
      read(nunitd,spc3,iostat=ios)
      if (ios /= 0) then
       write(*,*) 'could not read namelist spc3.in'
      else
       close(nunitd)
      endif
    
      write(*,*) 'iband=',iband
      read (string,'(I10)') dumint
      dumintminus = dumint - 1
      stringminus = trim(int2char(dumintminus))
      if (iband .eqv. .true.) then
       unit_prevend = trim(trim(prefix)//"/dat/"
     .               //"specyl_"//trim(stringminus)//"_end.dat")
      else
       unit_prevend = " "
      endif
      write(*,*) 'debug_files ',ibein, trim(unit_prevend)

!#! Marco, inizialize arrays
      call init_arrays
!#! Marco, inizialize variables
      call init_variab

      write(6,*) ' SIMULATION SETUP '
      write(6,*) ' eta0 = ', eta0,'  mue = ',mue
      write(*,*) ' aspect ratio',rr
      write(6,*) ' lundquist = ',1.d0/eta0,'  alet ',alet
      write(6,*) ' prandtl = ',mue/eta0,'  almu ',almu
      write(6,*) ' hart = ',1.d0/(eta0*sqrt(mue/eta0))
      write(6,*) ' itend = ',itend,' dt =',dt
      write(6,*) ' ntp= ',itend/ib
      write(6,*) ' zeit_f (tauA) =  ',itend*dt
      write(6,*) ' correzione flusso:  itor ',itor
      write(6,*) ' ipr = ', ipr, ' lx = ',lx
      write(6,*) ' numero di modi iz ',iz
      write(6,*) ' eps1 modi ',eps1
      write(6,*) ' ivac= ',ivac
      write(6,*) ' inopert= ',inopert
      write(6,*) ' inoconv= ',inoconv
      write(6,*) ' i2d= ',i2d
      write(6,*) ' ippcd= ',ippcd
      write(6,*) ' ilinear= ',ilinear
      write(6,*) ' lambda (di norma deve essere 1) = ',lambda
      write(6,*) ' c1 (mettere >0 per regione centrale calda) = ',c1
      write(6,*) ' GAET (di norma deve essere 1) = ',GAET
      write(6,*) ' irrotational MP? ',irrot
      write(6,*) ' MP = ',mp_perc,' MP phase = ',ph_mp
      write(6,*) ' rotating MP frequency = ',br_pert_freq
      write(6,*) ' Momentum source = ',m0_momsour
!      write(6,*) ' prefix = ',trim(prefix)
!      write(6,*) ' prefix = ',unit_all
      write(6,*) '   '
      call flush(6)

c    parti modificate dove compare: bza*   
c    introdotto - bzaini e bzafin (in cyl1*.for)
c               - scelta b.c.
c               - ecluso controllo sul flusso
   
      IF (I2D.EQ.1) THEN
        HELI2D = N2D*1./M2D
        nanf = mm * n2d
      ENDIF

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
       DT2 = DT*0.5d0
!#! Marco. Si fa un passo di ADT=ALF*DT (ALF=0.6) nella prima legge di faraday e poi
!un passo di DT2=0.5DT nella seconda legge di faraday. perche?
!#! Con Daniele il 12/09/2014 abbiamo capito che il 05*DT nella seconda
!legge di Faraday è solo una conseguenza della discretizzazione adottata
!(c'era un fattore due e lo hanno infilato nel DT), quindi ritorno a
!DT2=DT*0.5
!sostituisco DT2=(1-ALF)*DT
!       DT2 = DT * (1-ALF)
       ADT = ALF*DT                                                     
       DTMUE = DT*MUE                                                   
       DTKAP = DT*KAPPA                                                 
       GAM1 = GAMMA-1.0                                                 
       P0 = 47.0/(180.0*Q0**2)                                          
                                                                        
      DO 444 IM=0,MY                                                    
      M= MM(IM)                                                         
      DO 444 N=NANF(IM),NANF(IM)+NZ(IM)-1                               
                                                                        
      WRITE(6,'(a,i0,a,i0,a,i0,a,i0)') 'M= ',M,' N= ',N,
     .     ' JANZ(M,N)=',JANZ(M,N),' MY=',MY
                                                                        
444   CONTINUE                                                          
      call flush(6)

       DO 1012 IM=0,MY                                                    
       M = MM(IM)                                                       
!       write (6,*) 'm',IM,na1(m),ne(m)
      call flush(6)
        DO 1012 N=NA1(M),NE(M)    
!       DO 1012 N=NANF(IM),NANF(IM)+NZ(IM)-1                               
          nj=0
          jmn = JANZ(M,N)                                                  
          if (m.eq.0) jmn = -jmn ! serve altrimenti ho jmn negativi!!!
!          write (6,*) 'm,n,mnj',m,n,jmn
         DO 1056 MS=-MM(MY)+M,MM(MY)                                        
! non cambia nulla togliendo +M !!! (ma probabilmente fa più interazioni a vuoto)
!       DO 1056 MS=-MM(MY),MM(MY)                                        
          IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) cycle
          NS0 = NA(MS)
          NS1 = NE(MS)
          IF (I2D.EQ.1) THEN
             NS0 = MS*HELI2D
             NS1 = MS*HELI2D  
          ENDIF
          DO 1056 NS=NS0,NS1
          IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) cycle
          J = JANZ(MS,NS)                                                  
          J1 = JANZ(M-MS,N-NS)                                             
          CONVJ(jmn,nj)=j
          CONVJ1(jmn,nj)=j1
          CONVMS(jmn,nj)=ms
          CONVNS(jmn,nj)=ns
          nj=nj+1
!       write (6,112) ms,ns,m-ms,n-ns,j,j1,nj
 1056    CONTINUE
        CONVNJ(jmn)=nj
 1012   CONTINUE                                                         
 112    FORMAT('ms,ns,m-ms,n-ns,j,j1,nj',7I10)

!                                                                       
! ................ RESONANCE .....................                      
!                                                                       
!       IF (IPR.GE.0) THEN
! daniele, 9/11/2006
! modifica messa per evitare casini quando faccio simulazioni 1d (IZ=0)
! NB rimane RS = 0.4 ma non mi cambia nulla! (conta solo se è zero o non 0)
!       IF (IPR.GE.0.AND.IZ.GT.0) THEN
!         M = MM(1)                                                      
!         N = NANF(1)
!         RS = -1.                                                       
!         NRAD = 100                                                     
!         DRAD = 1./(NRAD - 1)                                           
!         RAD = 0.0                                                      
!         CALL EQUIL(RAD,BZ00,BT00,P010)                                 
!         F0 = M*THETA0 + N/RR*BZ00                                      
!         DO 785 IRAD=2,NRAD                                             
!              RAD = (IRAD - 1)*DRAD                                     
!              CALL EQUIL(RAD,BZ00,BT00,P010)                            
!              F01 = M/RAD*BT00 + N/RR*BZ00                              
!              DUM = F0*F01
!              IF (DUM.LE.0.0) THEN                                      
!                   RS = RAD - F01*DRAD/(F01 - F0)                       
!                   GO TO 786                                            
!              END IF                                                    
!              F0 = F01                                                  
! 785     CONTINUE                                                       
! 786     CONTINUE                                                       
!       END IF                                                           
                                                                        
! ....................MESH...............................               
                                                                        
! daniele, settembre 2006
! tolta perchè dava fastidio (nel caso tokamak)
! quanto q0>1
!       IF (IPR.EQ.-1) RS = SQRT(1.0/Q0-1.0)*C1+0.05                     
       dx = dble(1.d0/lx)
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
!       CALL RTNI( real(X(IX)),real(FUNC),real(DERF)
!     .             ,FCT,real(XST),1.0-8,40,real(IER1))   
!#! Marco, RTNI & company sono troppo troppo obsolete e non funzionano più
       x(ix) = dble(ix) * dx
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
!       WRITE(6,*)'ZDR: ',ZDR
!       WRITE(6,*)'X,X2,ZS2,ZSQ,ZZZ,ZZ1,Z2S2,ZS22,ZSQ2,ZZZ2',ZS,'ZS'    
!       DO 720 IX=1,100
!        WRITE(6,1003) X(IX),X2(IX),ZS2(IX),ZSQ(IX),ZZZ(IX),            
!     .                 ZZ1(IX),Z2S2(IX),ZSQ2(IX),ZZZ2(IX)               
!720    CONTINUE                                                        
!1003   FORMAT(10(f10.3))                                            
                                                                        
                                                                        
       WRITE (6,*) ' MESH FATTA '                                      

!#! Marco, compute momentum source
       call get_momsour(momsour)
 
!...................... INITIAL EQUILIBRIUM .....................      
                                                                        
!       if (ipr .eq. -2) then                               
!        call compute_equil(x,x2,ebt0,ebz0,evr0,bzwall,btwall)
!
!          write (*,*) 'ifuori  Bz,Bt(LX)=',ebz0(lx),ebt0(lx)
!       VRLX = evr0(LX)
!!#! Marco, assegna i campi calcolati dalla subroutine compute_equil
!       DO 435 IX=0,LX
!          BZ(IX,J0) = ebz0(IX)
!          BT(IX,J0) = ebt0(IX)
!          VR(IX,J0) = evr0(IX)
! 1002     FORMAT (1X,7E12.3)                                               
! 435   CONTINUE                                                          
!
!          write(*,*) 'outside equil', bzwall,btwall
!!          bzwall = ohmbz(LX)
!!          btwall = ohmbt(LX)
!!       ebz0(lx) = 2.d0*bzwall - ebz0(lx-1)
!!       ebt0(lx) = 2.d0*btwall - ebt0(lx-1)
!!       BZ(LX,J0) = ebz0(LX) 
!!       BT(LX,J0) = ebt0(LX) 
!!       VR(LX,J0) = evr0(LX)
!       BZ(-1,J0) = BZ(0,J0) 
!       BT(-1,J0) = -BT(0,J0)
!           
!!       A = DT * BETA * BZ(0,J0)                                       
!!       A22= A ** 2.                                                   
!
!       else
!        write(*,*) 'equilibrium with IPR=',ipr,'  is not implemented!'
!        stop
!       endif
!!#! Marco, bug
!
!!       BZ0(-1) =  BZ0(0)                                              
!!       BT0(-1) = -BT0(0)
!       do k=-1, lx
!        write(*,*) k,bz(k,0),bt(k,0),vr(k,0)
!       enddo
!       WRITE (6,*) ' BZ0(LX) = ',ebz0(LX),' BT0(LX) = ',ebt0(LX)
!       write(*,*) 'iband',IBAND
!       write(*,*) 'vrlx',VRLX
!       call flush(6)
!#!Marco, fine equilibrio nuovo

! Daniele, 09/05/2007
! implemento il calcolo dell'equilibrio screw pinch ohmico come in PIXIE3D
! qui setto alpha, se poi voglio modificare theta mi cambia B0..
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
!#!Marco, calcolo vz0 da momentum source in direzione z
!#! questa è analitica, per calcolare le altre componenti è da fare
!numericamente
          do ix=0,lx
             vz0(ix) = - momsour(ix,2) / (4.d0 * MUE) *
     .                   ( x2(ix)**2.d0 - 1.d0 )
!             write(*,*) 'mrd0 ',x2(ix), vz0(ix)
          enddo
          VRLX = VR0(LX)

!#!Marco, check on the core value of vz0
          if ( vz0(0) .gt. 0.5d0) then
           write(*,*) 'Momentum source too high. Velocity is greater
     .                 than vA/2! Stopping the code.'
           stop
          endif

          write (*,*) 
          write (*,*) 'Pinch equilibrium features:'
          write (*,*) '  Bz(r=0)=',ohmbz(0)
          write (*,*) '  Bth(r=a)=',ohmbt(LX)
          write (*,*) '  Bz avg=',bz_avg
          write (*,*) '  Theta  =',ohmbt(LX)/bz_avg
          
          write (*,*) '  I0     =',2.*pi*x(LX)*ohmbt(LX)
          write (*,*) '  E0/eta =',e0_over_eta
          write (*,*) '  E0 =',E0
          write (*,*) '  T_flux =',pi*x(LX)**2.*bz_avg
          write (*,*) '  Momsour =',m0_momsour
          write (*,*) '  Vz0(0) =',vz0(0)

!!          do ix=0,LX-1
!          do ix=0,LX
!!             write(*,1001) ix,x(ix),ohmbt(ix),ohmbz(ix)
!             write(*,1000) ix,x2(ix),BT0(ix),BZ0(ix),VR0(ix)
!          enddo
!!          write(*,1000) LX,x(LX),ohmbt(LX),ohmbz(LX),VR0(LX)
!          write (*,*) 'i,r,q0'
!          do ix=0,LX
!             write(*,1000) ix,x(ix),
!     .              (x(ix)*ohmbz(ix))/(RR*ohmbt(ix))
!          enddo
!          call flush(6)
!
! 1000     FORMAT(I3,' ',f6.3,' ',e12.5,' ',e12.5,' ',e12.5)
! 1001     FORMAT(I3,' ',f6.3,' ',e12.5,' ',e12.5)

!#! Marco, fine calcolo equilibrio ohmico
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
!             CALL EQUIL(X2(IX),BZ0(IX),BT0(IX),P01(IX))                    
! Daniele, 16/05/2007
!             VR0(IX) = 0.
          END IF
          BZ(IX,J0) = BZ0(IX)
          BT(IX,J0) = BT0(IX)
          VR(IX,J0) = VR0(IX)
!#! Marco, momentum source in z direction
          VZ(IX,J0) = vz0(IX)
!          WRITE (*,*) 'mrd ',X2(IX), vz(IX,J0)                     
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
!          CALL EQUIL(X(LX),BZWALL,BTWALL,P0WALL)
! Daniele, 16/05/2007
!          VR0(LX) = 0.
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
!!#!Marco, to consider a momentum source in z direction
       VZ(-1,J0) = vz0(0)
!       VZ(LX,J0) = - vz0(LX)
!         WRITE(6,1002) X(LX),BZ(LX,J0),BT(LX,J0)                       
       BZ(-1,J0) = BZ(0,J0)                                           
       BT(-1,J0) = -BT(0,J0)                                          
                                                                        
       A = DT * BETA * BZ(0,J0)                                       
       A22= A ** 2.                                                   

       WRITE (6,*) ' EQUILIBRIO INIZIALE FATTO: '
       WRITE (6,*) ' BZ0(LX) = ', bz0(LX),' BT0(LX) = ',bt0(LX)
       call flush(6)
                                                                       
!       WRITE (6,*) j0
!       do ix=-1,100 
!        write(*,*) X2(IX), BZ(IX,J0), BT(IX,J0), vr(ix,j0)
!       enddo
!#!Marco, fine equilibrio vecchio


!!...................... INITIAL  DENSITY     .....................      
                                                                        
!      DO 706 IM=0,MY                                                   
!      M = MM(IM)                                                       
!       DO 706 N=NA1(M),NE(M)                                           
!       J= JANZ(M,N)                                                    
!       J1=JANZ(-MM(IM),-N)                                             
!         DO 706 IX=0,LX                                                
!                                                                       
!         DE(IX,J) = (0.,0.)                                            
!         PR(IX,J) = (0.,0.)                                            
!         IF( J .EQ. 0 ) DE(IX,J)= (1.,0.)                              
!         DE(IX,J1) = CONJG( DE(IX,J) )                                 
!                                                                       
!         IF( J .EQ. 0) THEN                                            
!         PRRE= 2./3.                                                   
!         IF( IX .EQ. 0 ) PRR = 1.                                      
!         IF(IX.GT.0) PRR = (1.+(ALET-1.)*X(IX)**BEET)**PRR             
!         PRR = 1./PRR                                                  
!         PR(IX,J)= (1.,0.)                                             
!         PR(IX,J)= PRR * PR(IX,J) * 0.5                                
!         WRITE(6,*)'  ',IX,'   ', PRR                                  
!         ENDIF                                                         
!                                                                       
!         PR(IX,J1) = CONJG( PR(IX,J) )                                 
                                                                        
!         IF( IX .GE. 1)THEN                                            
!         AG = 3.14 * X(IX)                                             
!         AG2= 3.14 * X2(IX)                                            
!                                                                       
!         FIX  = 0.1 *(SIN(AG))**2.                                     
!         DFIX = 0.1 * PI * SIN(2.*AG2)                                 
!                                                                       
!         FIX  = 1.0                                                    
!         DFIX = 1.0 * 1.                                               
!                                                                       
!         FIX  = 0.1 * SIN(AG)                                          
!         DFIX = 0.1 * PI * COS( AG2)                                   
!         VRP(IX,J) = 0.                                                
!         VTP(IX,J) = 0.                                                
!         IF(M .EQ. 1 .AND. N .EQ. 0)THEN                               
!         VRP(IX,J) = FIX/(X(IX)*2.)                                    
!         VRP(IX,J) = FIX/ 2.                                           
!         VTP(IX,J) = I * DFIX / 2.                                     
!         IF( IX.EQ.0) WRITE(6,*) VTP(IX,J),J,'QUI'                     
!         ENDIF                                                         
!         ENDIF                                                         
!         VRP(IX,J1)= CONJG( VRP(IX,J) )                                
!         VTP(IX,J1)= CONJG( VTP(IX,J) )                                
!                                                                       
!706    CONTINUE                                                        
!...........................................................            
                                                                        
! - PROFILES OF THE DISSIPATION COEFFICIENTS                              
                                                                        
       call dissipation_profiles(mu1,mu2,etaa,eta2,kap)
       call flush(6)


       ZEIT = 0.0                                                       
       IBIB = 0                                                         
								
!#! Marco, read from disk only if IBAND .eq. TRUE
       IF(.NOT.IBAND) GO TO 88                                          
								
!..................... READ FROM DISK OR FRONT-END ...............      
								
       WRITE(6,*)'INIZIO LETTURA'
       								
!       REWIND IBEIN                                                     
       write(*,*) trim(unit_prevend),ibein
       open(unit=ibein,file=trim(unit_prevend),action='read'
     .        ,form='unformatted')

       WRITE(6,*)'ZEIT=',ZEIT,'   LY',LY,'    IZZ',IZZ
       IF (IBAND1 .EQ. 0) then
        READ(IBEIN) ZEIT,LY,IZZ
!        READ(IBEIN) ZEIT
!        READ(IBEIN) LY,IZZ
       endif
       WRITE(6,*)'ZEIT=',ZEIT,'   LY',LY,'    IZZ',IZZ
       call flush(6)
       								
       IF( IZZ .NE. IZ ) then
       WRITE(6,*) 'NUMERO DIVERSO DI MODI ELABORATI:IZZ/IZ'
       endif
       call flush(6)
       
       DO 250 J=0,IZZ
       IF(IBAND1.EQ.0) then
        READ(IBEIN) VT(-1,j), VZ(-1,j), BT(-1,j), BZ(-1,j)
       endif
         VT(-1,-J) = CONJG(VT(-1,J))                                   
         VZ(-1,-J) = CONJG(VZ(-1,J))                                   
         BT(-1,-J) = CONJG(BT(-1,J))                                   
         BZ(-1,-J) = CONJG(BZ(-1,J))                                   
       DO 250 IX=0,LY                                                   
         IF(IBAND1.EQ.0) then
          READ(IBEIN) VR(ix,j), VT(ix,j), VZ(ix,j),
     .                 BR(ix,j), BT(ix,j), BZ(ix,j) 
         endif
         VR(IX,-J) = CONJG(VR(IX,J))                                   
         VT(IX,-J) = CONJG(VT(IX,J))                                   
         VZ(IX,-J) = CONJG(VZ(IX,J))                                   
         BR(IX,-J) = CONJG(BR(IX,J))                                   
         BT(IX,-J) = CONJG(BT(IX,J))                                   
         BZ(IX,-J) = CONJG(BZ(IX,J))                                   
       
      
250   CONTINUE          
        close(ibein)
       
       !   THE PERTURBATIONS (J=-1,1) ARE NORMALIZED TO EPS1                   
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
       								
       !#! Marco, if starting from initial equilibrium begin from here
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
       WRITE(6,*) 'x= ',X2(50),' bz ',BZ(50,0),' bt ',BT(50,0)
       WRITE(6,*) 'x= ',X2(50),' br ',BR(50,0)
       
       WRITE(6,*) '   REAL(BZ(LX-1,J0)) =',REAL(BZ(LX-1,J0)),
     .        '   REAL(BZ(LX,J0)) = ', REAL(BZ(LX,J0))
       WRITE(6,*) '   REAL(BT(LX-1,J0)) =',REAL(BT(LX-1,J0)),
     .        '   REAL(BT(LX,J0)) = ', REAL(BT(LX,J0))
       !       WRITE(6,*) '   0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0))) = ',
       !     .        0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0)))
       WRITE(6,*) '   TH =',YTH, '   F = ', YFITP,'   zeit = ',zeit
       WRITE(6,*) '   flusso toroidale= ', TORF1
       WRITE(6,*) '   Qo= ', YQITP
       
       call flush(6)
								

!                    perturbazione per il run :
!#! Marco, 17/02/2017
       write(*,*) 'inopert=',inopert
       if (inopert .eq. 0) then
        call initial_pert
       endif

       call flush(6)
       write(*,*) 

!#! Marco, 20/02/2017
!#! subroutine che impone le condizioni al contorno del campo magnetico
!all'inizio di ogni tranche di simulazione
!#! La subroutine tiene in conto i seguenti casi: MP su più modi, MP con
!fase iniziale scelta dall'utente, MP rotanti.

       call apply_mp_init


! daniele, 10/12/2008
! scrivo anche l'istante t=0 che contiene l'equilibrio più la perturbazione assegnata
! oppure se stavo partendo da un'altra simulazione, ho la configurazione da cui parto
!#! Marco 06/08/2014 I change how data are written on disk. The whole
!complex number will be written, not just its real or imaginary part.
!........................................ WRITE ON DISK                 

!#! Marco, apri il file specyl_?_all.dat
       open(unit=ibout1,file=trim(unit_all),action='write'
     .        ,form='unformatted')
!#! Marco, apri il file specyl_?_end.dat        
        open(unit=ibaus,file=trim(unit_end),action='write'
     .        ,form='unformatted')

       write(*,*) 'inizio scrittura',zeit,lx,iz
       WRITE(IBOUT1)  ZEIT,LX,IZ
!       WRITE(IBOUT1)  ZEIT
!       WRITE(IBOUT1)  int(LX),int(IZ)
       DO 1252 J=0,IZ                                                    
        WRITE(IBOUT1) VT(-1,j), VZ(-1,j), BT(-1,j), BZ(-1,j)
       
       DO 1252 IX=0,LX
         WRITE(IBOUT1) VR(ix,j), VT(ix,j), VZ(ix,j)
     .                ,BR(ix,j), BT(ix,j), BZ(ix,j)
1252   CONTINUE                                                         
       call flush(IBOUT1)
       call flush(6)
       
       !.................. MAIN LOOP STARTS HERE...............................
       								
       WRITE(6,*)'INIZIO MAIN LOOP'                                    

c..... bza definizione:
       bzaini = BZ(LX,0)
!             if (ippcd.ne.0) then
!!       bzafin = bzaini
!!                write(6,*) ' ippcd = ',ippcd
!                write(6,*) ' bza_ini = ',bzaini
!                write(6,*) ' bza_fin = ',bzafin
!                call flush(6)                
!             endif
								
!             write(*,*) 'va',j_mp
       DO 50 IT=1,ITEND                                                 
         ZEIT = ZEIT+DT                                                   
       
       acontr = (bt(50,j0) * conjg(bt(50,j0)))
       afluct = (bt(50,1) * conjg(bt(50,1)))
       if(acontr.gt.1.e5) write(6,*)'acontr fuori contr ! it=',it
     .                         ,acontr
       if(acontr.gt.1.e5) go to 5000
       
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
       								
       J = JANZ(MPLO(2),NPLO(2))                                        
       VRALT = VR(IG,J)                                                 
       VZALT = VZ(IG,J)                                                 
       								
!............... B UND P - HALBSCHRITT (mezzo step) .................................
       !#! Marco, prima legge di Faraday
       								
       DO 12 IM=0,MY                                                    
       M = MM(IM)                                                       
       DO 12 N=NA1(M),NE(M)                                             
         SUM0=(0.0d0,0.0d0)
         SUM1=(0.0d0,0.0d0)
         SUM2=(0.0d0,0.0d0)
         SUM3=(0.0d0,0.0d0)
         SUM8=(0.0d0,0.0d0)
       
       nj = 0
       jmn = janz(m,n)
       if (m.eq.0) jmn = -jmn ! serve altrimenti ho jmn negativi!!!
       !#! Marco, debug
       !       if (m .eq. 0 .and. n .eq. 0 .and. lx .eq. 50) then
       !         write(*,*) 'it=',it
       !        write(*,'(a,6i4)') 'm,n,jmn,convnj',m,n,jmn,convnj(jmn)
       !       endif
!$OMP PARALLEL DO DEFAULT(SHARED)
!#! Marco
!$OMP& PRIVATE(IAM,js,IX,J,J1,m,n,j_mp)
!$OMP& SCHEDULE(STATIC)
!$OMP& REDUCTION(+:SUM0,SUM1,SUM2,SUM3,SUM8)
       DO 56 js=0,convnj(jmn)-1
       IAM=OMP_GET_THREAD_NUM()
       j = convj(jmn,js)
       j1 = convj1(jmn,js)
       !#! Marco, debug
       !       if (m .eq. 0 .and. n .eq. 0 .and. lx .eq. 50) then
       !        write(*,'(a,6i4)') 'm,n,jmn',m,n,jmn,js,j,j1
       !       endif
       !       write (6,111) iam,ms,ns,m-ms,n-ns
       !       write (6,111) iam,jmn,js,j,j1
       DO 58 IX=0,LX-1                                                  
       SUM3(IX) = SUM3(IX)+VT(IX,J)*BZ(IX,J1)-VZ(IX,J)*BT(IX,J1)        
58     CONTINUE                                                          
       DO 26 IX=0,LX                                                    
       SUM0(IX) = SUM0(IX)+BR(IX,J)*(VT(IX-1,J1)+VT(IX,J1))             
       SUM1(IX) = SUM1(IX)+VR(IX,J)*(BT(IX-1,J1)+BT(IX,J1))             
       SUM2(IX) = SUM2(IX)+BR(IX,J)*(VZ(IX-1,J1)+VZ(IX,J1))             
       SUM8(IX) = SUM8(IX)+VR(IX,J)*(BZ(IX-1,J1)+BZ(IX,J1))             
       !#! Marco, debug
       !       if (m .eq. 0 .and. n .eq. 0 .and. ix .eq. 50) then
       !        write(*,*)  'an= ',js,j,j1,j_mp
       !!     .           ,sum0(ix), sum1(ix),vr(ix,j)
       !!     .           ,br(ix,j), vr(ix,j)
       !!     .           ,vt(ix,j1), bt(ix,j1)
       !     .           ,br(ix,j_mp),br(ix,1)
       !       endif
26     CONTINUE                                                          
!   57 CONTINUE                                                          
56     CONTINUE
!$OMP END PARALLEL DO
       
       !       write (6,*) m,n,jmn,convnj(jmn)-1,js
       !  111 FORMAT('iam,ms,ns,m-ms,n-ns',5I10)
       !  111 FORMAT('iam,jmn,js,j,j1',5I10)
       ! ACHTUNG MANCHE SUMMEN (ABL.) AUCH AM RAND NOTWENDIG                   
       DO 59 IX=0,LX-1                                                  
       DSUM2(IX)=(SUM1(IX+1)-SUM1(IX)-SUM0(IX+1)+SUM0(IX))
     .          *ZS22(IX)*0.5d0
       !       if (m .eq. 0 .and. n .eq. 0 .and. ix .eq. 50) then
       !         write(*,*) 'dsum2',dsum2(ix)
       !       endif
       
       DSUM3(IX) =((SUM2(IX+1)-SUM8(IX+1))*X(IX+1)-(SUM2(IX)-SUM8(IX))  
     .          *X(IX))*ZS22(IX)/X2(IX)*0.5d0
       !#! Marco
       !       if (m .eq. 1 .and. ix .eq. 50) then
       !        write(*,*) 'guardo dsum2,prima faraday ',zs22(ix),
       !     .          dsum2(ix),dsum3(ix)
       !!     .           sum1(ix),sum0(ix),dsum2(ix)
       !       endif
59     CONTINUE                                                          
       !#! Marco, perchè chiamano J1?????
       !       if (j .eq. j0) then
       !        write(*,*)  'primaprima solenoida,it= ',it,
       !     .             'btn=',btn(50,j0),' bt=',bt(50,j0)
       !       endif
       J = JANZ(MM(IM),N)                                               
       J1 = JANZ(-MM(IM),-N)                                            
       CFAK1 = N/RR*I                                                   
       DO 14 IX=0,LX-1                                                  
       CFAK = M/X2(IX)*I                                                
       BTN(IX,J) = BT(IX,J)+ADT*(+CFAK1*SUM3(IX)-DSUM2(IX))             
       BZN(IX,J) = BZ(IX,J)+ADT*(+DSUM3(IX)-CFAK*SUM3(IX))
14     CONTINUE                                                          
1112   format ('a,i0,a,e8.3,a,e8.3')
       !#!Marco
       !       if (j .eq. j0) then
       !        write(*,*)  'primaima solenoida2,it= ',it,
       !     .             'btn=',btn(50,j0),' bt=',bt(50,j0)
       !       endif
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
       
       !     ............B.C. PER BZ ALLA PARETE:                              
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
       !        BZN(LX,J) = BZ0(LX)                                            
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
       
       !    ........                                                           
66       IF(M .NE. 0) GO TO 68                                          
         BZN(-1,J) = BZN(0,J)                                           
68       J1 = JANZ(-M,-N)                                                 
       								
       
!  BR AUS BT +BZ ........................                               
       BRN(0,J) =(0.0,0.0)                                              
!#! Marco, calcolo anche il punto LX
!       DO 149 IX=1,LX-1                                                 
       DO 149 IX=1,LX
       BRN(IX,J)=(I/ZS22(IX-1)*(-BTN(IX-1,J)*M-BZN(IX-1,J)*N*X2(IX-1)   
     .           /RR)+BRN(IX-1,J)*X(IX-1))/X(IX)                                  
149     CONTINUE                                                          
       IF(M.NE.1) GO TO 67                                               
       !      BTN(0,J) = -BRN(0,J)                                             
       BRN(0,J) = (4.0*BRN(1,J)-BRN(2,J))*0.333333333                   
       BTN(-1,J) = -BTN(0,J)+2.0*I*BRN(0,J)                             
67     IF(J.EQ.0) GO TO 12                                              
       DO 98 IX=0,LX                                                    
       BRN(IX,J1) = CONJG(BRN(IX,J))                                    
98     CONTINUE                                                          
       DO 401 IX=-1,LX                                                  
       BTN(IX,J1) = CONJG(BTN(IX,J))                                    
       BZN(IX,J1) = CONJG(BZN(IX,J))                                    
401     CONTINUE            
12     CONTINUE                                                          
       !#! Marco
!        write(*,*),'solenoidalita1 ,it= ',it
!     .             ,'btn0=',btn(50,j0),' bt0=',bt(50,j0)
!     .             ,'brn0=',brn(50,j0),' br0=',br(50,j0)
!     .             ,'vtn0=',vtn(50,j0),' vt0=',vt(50,j0)
!     .             ,'brn=',brn(100,74),' br=',br(100,74)
!     .             ,'bzn=',bzn(100,74),' bz=',bz(100,74)
       								
       ! .................... J-BESTIMMUNG ................................... 
       								
       DO 410 IM=0,MY                                                   
       M = MM(IM)                                                       
       DO 410 N=NA1(M),NE(M)                                            
       FAK1 = N/RR                                                      
       J = JANZ(M,N)                                                    
       J1 = JANZ(-M,-N)                                                 
       DO 408 IX=0,LX-1                                                 
       JR(IX,J) = I*(M/X2(IX)*BZN(IX,J)-FAK1*BTN(IX,J))                 
408     CONTINUE                                                          
       !#! Marco, 02/10/2014, modifica fatta con Daniele per calcolare la
       !corrente anche a LX
       !       DO 409 IX=1,LX-1                                                 
       DO 409 IX=1,LX
       JT(IX,J) = I*FAK1*BRN(IX,J)-2.0*ZS2(IX)*(BZN(IX,J)-BZN(IX-1,J))  
       JZ(IX,J) = 2.0*ZS2(IX)/X(IX)*(X2(IX)*BTN(IX,J)-X2(IX-1)          
     .            *BTN(IX-1,J))-I*M/X(IX)*BRN(IX,J)                                
409     CONTINUE                                                          
       !       JT(LX,J) = (0.0,0.0)                                             
       !       JZ(LX,J) = (0.0,0.0)                                             
       JT(0,J) =(0.0,0.0)                                               
       JZ(0,J) =(0.0,0.0)                                               
       IF(M.NE.0) GO TO 434                                             
       !       JT(LX,J) = -2.*ZS2(LX)*(BZN(LX,J) - BZN(LX-1,J))                 
       !       JZ(LX,J) =  2.*ZS2(LX)*((2. - X2(LX-1))*BTN(LX,J) -              
       !     &                                          X2(LX-1)*BTN(LX-1,J))   
       JZ(0,J) = 2.*BTN(0,J)/X2(0)                                     
434    IF(M.NE.1) GO TO 420                                              
       JT(0,J) = JT(1,J)                                                
420    IF(J.EQ.0) GO TO 410                                              
       DO 437 IX=0,LX                                                   
       JR(IX,J1) = CONJG(JR(IX,J))                                      
437     CONTINUE                                                          
       DO 421 IX=0,LX                                                   
       JT(IX,J1) = CONJG(JT(IX,J))                                      
       JZ(IX,J1) = CONJG(JZ(IX,J))                                      
421     CONTINUE                                                          
410     CONTINUE                                                          
       								
!          write(*,*)  'medio= ',bz(90,j0)
! ............................DENSITY ..................................
								
!     DO 700 IM=0,MY                                                    
!      M = MM(IM)                                                       
!      DO 700 N=NA1(M),NE(M)                                            
!                                                                       
!       DO 701 IX = 0,LX                                                
!       SUM0(IX) = (0.0,0.0)                                            
!       SUM1(IX) = (0.0,0.0)                                            
!701     CONTINUE                                                       
!                                                                       
!       DO 702 MS = -MM(MY)+M,MM(MY)                                    
!      IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) GO TO 702                 
!        DO 703 NS = NA(MS),NE(MS)                                      
!      IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) GO TO 703               
!        J = JANZ(MS,NS)                                                
!        J1= JANZ(M-MS,N-NS)                                            
!          DO 704 IX=1,LX-1                                             
!                                                                       
!          VTIXJ= (VT(IX-1,J)+VT(IX,J))/2.                              
!          VZIXJ= (VZ(IX-1,J)+VZ(IX,J))/2.                              
!                                                                       
!     DENSITY :                                                         
!          S0 = - ZDR * 0.5 / X(IX)                                     
!          S1 = X(IX+1)*DE(IX+1,J1)*VR (IX+1,J)                         
!          S2 = X(IX-1)*DE(IX-1,J1)*VR (IX-1,J)                         
!          S3 = DE(IX,J1)                                               
!          S4 = I * M * VTIXJ  / X(IX)                                  
!          S5 = I * N * VZIXJ   / RR                                    
!          SUM0(IX) = S0 *(S1-S2) - S3 *(S4+S5) + SUM0(IX)              
!                                                                       
!    PRESSURE :                                                         
!          S0 = - ZDR * 0.5 / X(IX)                                     
!          S1 = X(IX+1)*PR(IX+1,J1)*VR (IX+1,J)                         
!          S2 = X(IX-1)*PR(IX-1,J1)*VR (IX-1,J)                         
!          S3 = PR(IX,J1)                                               
!          S4 = I * M * VTIXJ  / X(IX)                                  
!          S5 = I * N * VZIXJ   / RR                                    
!          S6 = -2./3. * PR(IX,J1) * (                                  
!    .    1./X(IX) * (X(IX+1)*VR(IX+1,J)-X(IX-1)*VR(IX-1,J)) *0.5 *ZDR  
!    .       + I * MS /X(IX) * VTIXJ + I * NS /RR * VZIXJ    )          
!          SUM1(IX) = S0 *(S1-S2) - S3 *(S4+S5) + S6 + SUM1(IX)         
!                                                                       
!704       CONTINUE                                                     
!703       CONTINUE                                                     
!702       CONTINUE                                                     
!                                                                       
!         J = JANZ(MM(IM),N)                                            
!         J1= JANZ(-MM(IM),-N)                                          
!         DO 705 IX=1,LX-1                                              
!                                                                       
!          AMX = 1.0                                                    
!          DENN(IX,J)=AMX*(DE(IX+1,J)+DE(IX-1,J))/2.+ (1.-AMX)*DE(IX,J) 
!    .                + ADT * SUM0(IX)                                  
!                                                                       
!          PRNN(IX,J)=AMX*(PR(IX+1,J)+PR(IX-1,J))/2.+ (1.-AMX)*PR(IX,J) 
!    .                + ADT * SUM1(IX)                                  
!                                                                       
!          DENN(IX,J1)= CONJG(DENN(IX,J))                               
!          PRNN(IX,J1)= CONJG(PRNN(IX,J))                               
!705       CONTINUE                                                     
!          IF ( MM(IM) .EQ. 0 ) THEN                                    
!          DENN(0,J)  = DENN(1,J)                                       
!          DENN(0,J1) = CONJG(DENN(0,J))                                
!          PRNN(0,J)  = PRNN(1,J)                                       
!         PRNN(0,J1) = CONJG(PRNN(0,J))                                 
!         ELSE                                                          
!         DENN(0,J)  = 0.0                                              
!         DENN(0,J1) = CONJG(DENN(0,J))                                 
!         PRNN(0,J)  = 0.0                                              
!         PRNN(0,J1) = CONJG(PRNN(0,J))                                 
!         ENDIF                                                         
!         IF( J .EQ. 0 )THEN                                            
!         DENN(LX,J) = DENN(LX-1,J)                                     
!         ELSE                                                          
!         DENN(LX,J) = 0.0                                              
!         ENDIF                                                         
!         DENN(LX,J1)=CONJG(DENN(LX,J))                                 
!         PRNN(LX,J) = PRNN(LX-1,J)                                     
!         PRNN(LX,J) = PR(LX,J)                                         
!         PRNN(LX,J1)=CONJG(PRNN(LX,J))                                 
!700       CONTINUE                                                     

!                                                                       
! ............  GANZSCHRITTSEMI-IMPLIZ. VERFAHREN ......................
!#! Marco, passo seminplicito
								
        DO 15 IM=0,MY                                                    
        M = MM(IM)                                                       
        DO 15 N=NA1(M),NE(M)                                             
        ampv = 1.e0
         FAK1 = N/RR 
         SUM0=(0.0d0,0.0d0)
         SUM1=(0.0d0,0.0d0)
         SUM2=(0.0d0,0.0d0)
         SUM3=(0.0d0,0.0d0)
         SUM4=(0.0d0,0.0d0)
         SUM5=(0.0d0,0.0d0)
         SUM6=(0.0d0,0.0d0)
         SUM7=(0.0d0,0.0d0)
         SUM8=(0.0d0,0.0d0)
         SUM9=(0.0d0,0.0d0)
         SUM10=(0.0d0,0.0d0)
         SUM11=(0.0d0,0.0d0)
         SUM12=(0.0d0,0.0d0)
         SUM13=(0.0d0,0.0d0)
         SUM14=(0.0d0,0.0d0)
         SUM15=(0.0d0,0.0d0)
         SUM16=(0.0d0,0.0d0)
         DSUM1=(0.0d0,0.0d0)
         DSUM2=(0.0d0,0.0d0)
         DSUM3=(0.0d0,0.0d0)
        
        jmn = janz(m,n)
        if (m.eq.0) jmn = -jmn ! serve altrimenti ho jmn negativi!!!
					     
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(IAM,js,IX,J,J1,ms,ns)
!$OMP& SCHEDULE(STATIC)
!$OMP& REDUCTION(+:SUM0,SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8,SUM9)
!$OMP& REDUCTION(+:SUM10,SUM11,SUM12,SUM13,SUM14,SUM15,SUM16)
        DO 17 js=0,convnj(jmn)-1
        IAM=OMP_GET_THREAD_NUM()
        j = convj(jmn,js)
        j1 = convj1(jmn,js)
        ms = convms(jmn,js)
        ns = convns(jmn,js)
        !       F1 = NS/RR                                                       
        DO 19 IX=1,LX-1                                                  
        SUM1(IX) = SUM1(IX)+(VR(IX+1,J)-VR(IX-1,J))*VR(IX,J1)            
        SUM2(IX)=SUM2(IX)+(VT(IX,J)+VT(IX-1,J))*(VT(IX,J1)+VT(IX-1,J1))  
        SUM3(IX) = SUM3(IX)+MS*VR(IX,J)*(VT(IX,J1)+VT(IX-1,J1))          
        SUM4(IX) = SUM4(IX)+NS*VR(IX,J)*(VZ(IX,J1)+VZ(IX-1,J1))          
        SUM12(IX) = SUM12(IX)+JT(IX,J)*(BZN(IX,J1)+BZN(IX-1,J1))         
        SUM0(IX) = SUM0(IX)+JZ(IX,J)*(BTN(IX,J1)+BTN(IX-1,J1))           
19     CONTINUE                                                          
        DO 27 IX=1,LX-2                                                  
        SUM5(IX)=SUM5(IX)+(VT(IX+1,J)-VT(IX-1,J))*
     .          (VR(IX,J1)+VR(IX+1,J1))
        SUM6(IX)=SUM6(IX)+VT(IX,J)*(VR(IX,J1)+VR(IX+1,J1))               
        SUM7(IX)=SUM7(IX)+MS*VT(IX,J)*VT(IX,J1)                          
        SUM8(IX)=SUM8(IX)+NS*VT(IX,J)*VZ(IX,J1)                          
        SUM9(IX)=SUM9(IX)+(VZ(IX+1,J)-VZ(IX-1,J))*
     .          (VR(IX,J1)+VR(IX+1,J1))
        SUM10(IX)=SUM10(IX)+MS*VZ(IX,J)*VT(IX,J1)                        
        SUM11(IX)=SUM11(IX)+NS*VZ(IX,J)*VZ(IX,J1)                        
27     CONTINUE                                                          
        DO 423 IX=0,LX-1                                                 
        SUM13(IX) = SUM13(IX)+(JZ(IX+1,J)+JZ(IX,J))*(BRN(IX+1,J1)
     .  + BRN(IX,J1))                                                      
        SUM14(IX) = SUM14(IX)+JR(IX,J)*BZN(IX,J1)                        
        SUM15(IX) = SUM15(IX)+JR(IX,J)*BTN(IX,J1)                        
        SUM16(IX) = SUM16(IX)+(JT(IX+1,J)+JT(IX,J))*(BRN(IX+1,J1)
     .  + BRN(IX,J1))                                                      
423     CONTINUE                                                          
        !   18     CONTINUE                                                          
17     CONTINUE                                                          
!$OMP END PARALLEL DO
        ! daniele, 24/11/2006
        ! aggiungo per escludere il termine convettivo
        IF (INOCONV.EQ.1) GO TO 399

        DO 429 IX=1,LX-1                                                 
        DSUM1(IX) = ZS2(IX)*SUM1(IX)-0.25/X(IX)*SUM2(IX)+
     .              0.5*I*(SUM3(IX)/X(IX)+SUM4(IX)/RR)                                              
429     CONTINUE                                                          
        
        DO 424 IX=1,LX-2                                                 
        DSUM2(IX)=(SUM5(IX)*ZS22(IX)*0.5+SUM6(IX)/X2(IX))*0.5
     .          +I*(SUM7(IX)/X2(IX)+SUM8(IX)/RR)                                             
        DSUM3(IX) = SUM9(IX)*ZS22(IX)*0.25+I*(SUM10(IX)
     .  /X2(IX)+SUM11(IX)/RR)                                             
424     CONTINUE                                                          
        
399     DSUM2(LX-1) = DSUM2(LX-2)                                        
        DSUM3(LX-1) = DSUM3(LX-2)                                        
        DSUM2(0) =(0.0d0,0.0d0)
        DSUM3(0) =(0.0d0,0.0d0)                                              
        IF(M.NE.1) GO TO 411                                             
        DSUM2(0) = DSUM2(1)                                              
411     IF(M.NE.0) GO TO 412                                             
        DSUM3(0) = DSUM3(1)                                              
412     J = JANZ(M,N)                                                    
        J1 = JANZ(-M,-N)                                                 
        								
        DO 20 IX=1,LX-1                                                  
        !        RHOX = REAL (DENN(IX,0))                                       
        !        IF( IACC .NE. 1 ) RHOX = 1.                                    
         KMN1 = (M/X(IX))**2 + (N/RR)**2                                
        !#! Marco, radial component of the rhs of equation of motion
         R1(IX) = VR(IX,J)                                              
     .   - A22*( ZDR**2*(VR(IX+1,J)-2.*VR(IX,J)+VR(IX-1,J))             
     .        + ZDR*0.5/X(IX)*(VR(IX+1,J)-VR(IX-1,J))                   
     .        - (KMN1+1./ X(IX)**2)*VR(IX,J)                            
     .        - I*M/X(IX)**2*(VT(IX,J)+VT(IX-1,J)) )                    
     .   - DT*( DSUM1(IX) +  0.5*(-SUM12(IX) + SUM0(IX)) )              
        !#! Marco, axisymmetric momentum source
         if (m .eq. 0 .and. n .eq. 0) then
          r1(ix) = r1(ix) + dt * momsour(ix,0) 
         endif
        !    5   - DT*( DSUM1(IX) +  0.5*(-SUM12(IX) + SUM0(IX)                 
        !    6          + BPRES * (PRNN(IX+1,J)-PRNN(IX-1,J)) * ZDR             
        !    7           )/RHOX )                                               
20     CONTINUE                                                          
        								
        DO 430 IX=0,LX-1                                                 
        !        RHOX2 = ( REAL(DENN(IX,0))+REAL(DENN(IX+1,0)) )/2.             
        !        IF ( IACC .NE. 1 ) RHOX2 = 1.                                  
         KMN2 = (M/X2(IX))**2 + (N/RR)**2                               
        !#! Marco, poloidal component of the rhs of equation of motion
         R2(IX) = VT(IX,J)                                              
     .   - A22 * ( ZDR**2/X2(IX)*                                       
     .     (X(IX+1)*(VT(IX+1,J)-VT(IX,J))-X(IX)*(VT(IX,J)-VT(IX-1,J)))  
     .     - (KMN2+1./X2(IX)**2)*VT(IX,J)                               
     .     + I*M/X2(IX)**2*(VR(IX+1,J)+VR(IX,J))  )                     
     .  - DT * ( DSUM2(IX)  -SUM13(IX)*0.25+SUM14(IX)  )                
        !#! Marco, axisymmetric momentum source
         if (m .eq. 0 .and. n .eq. 0) then
          r2(ix) = r2(ix) + dt * momsour(ix,1) 
         endif
        !    5  - DT * ( DSUM2(IX) +( -SUM13(IX)*0.25+SUM14(IX)                 
        !    6    + BPRES * I * M /X2(IX) *(PRNN(IX-1,J)+PRNN(IX+1,J))/2.       
        !    7           )/RHOX2 )                                              
        !#! Marco, axial component of the rhs of equation of motion
         R3(IX) = VZ(IX,J)                                              
     .   - A22 * ( ZDR**2/X2(IX)*                                       
     .      (X(IX+1)*(VZ(IX+1,J)-VZ(IX,J))-X(IX)*(VZ(IX,J)-VZ(IX-1,J))) 
     .      - KMN2*VZ(IX,J)  )                                          
     .   - DT * ( DSUM3(IX)  -SUM15(IX)+SUM16(IX)*0.25 )                
        !#! Marco, axysimmetric momentum source
         if (m .eq. 0 .and. n .eq. 0) then
          r3(ix) = r3(ix) + dt * momsour(ix,2) 
         endif
        !    4   - DT * ( DSUM3(IX) +( -SUM15(IX)+SUM16(IX)*0.25                
        !    5      + BPRES * I * N / RR *(PRNN(IX-1,J)+PRNN(IX+1,J))/2.        
        !    6           )/RHOX2 )                                              
430     CONTINUE                                                          
        								
         DO 801 L=1,2                                                   
              FF(L,0) = (0.0d0,0.0d0)                                       
         DO 801 K=1,2                                                   
              EE(L,K,0) = (0.0d0,0.0d0)                                     
801     CONTINUE                                                       
        								
         A12 = ampv * MU2(0)*DT + A22                                          
        !        RHO02 = (REAL(DENN(0,0))+REAL(DENN(1,0)))/2.                   
        !        IF( IACC .NE. 1 ) RHO02 = 1.                                   
        !        A12 = MU2(0)*DT/RHO02 + A22                                    
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
        !         RHOX1= REAL(DENN(IX,0))                                   
        !         RHOX2= (REAL(DENN(IX,0))+REAL(DENN(IX+1,0)))/2.           
        !         IF ( IACC .NE. 1 ) RHOX1= 1.                              
        !         IF ( IACC .NE. 1 ) RHOX2= 1.                              
          A11 = ampv * MU1(IX)*DT + A22                                    
          A12 = ampv * MU2(IX)*DT + A22                                    
        !         A11 = MU1(IX)*DT/RHOX1 + A22                              
        !         A12 = MU2(IX)*DT/RHOX2 + A22                              
                                                                            
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
        !#! Marco
        !         if (m .eq. 0 .and. n .eq. 0 .and. ix .eq. 50
        !     .      .and. j .eq. 0 .and. j1 .eq .0) then
        !          write(*,*) 'semimpl0', cc(:,:,ix),dd(:,ix)
        !         endif
                                                                                
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
        !#! Marco
        !         if (m .eq. 0 .and. n .eq. 0 .and. ix .eq. 50
        !     .      .and. j .eq. 0 .and. j1 .eq .0) then
        !          write(*,*) 'semimpl00', av3(1), av3(2),
        !     .                ee(2,2,ix)
        !         endif
             AAA(IX) = -ALJ(IX)/(BETJ(IX) + GAMJ(IX)*AAA(IX-1))        
             BBB(IX) = (R3(IX) - GAMJ(IX)*BBB(IX-1))/                  
     .                 (BETJ(IX) + GAMJ(IX)*AAA(IX-1))        
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
     .       /(1. + EE(2,2,LX-1))
             UU(2,LX) = - UU(2,LX-1)
             UU(1,LX-1) = EE(1,1,LX-1)*UU(1,LX) + EE(1,2,LX-1)*UU(2,LX)
     .       + FF(1,LX-1)
        
             VZN(LX-1,J) = BBB(LX-1)/(1. + AAA(LX-1))
             VZN(LX,J) = - VZN(LX-1,J)
        
        !         DO 807 IX=LX-1,0,-1
        ! daniele, 02/12/2008
        ! modifiche alle B.C. di vtheta e vz, per averle 0. esattamente a r=1.
             DO 807 IX=LX-2,0,-1
               UU(1,IX) = EE(1,1,IX)*UU(1,IX+1) + EE(1,2,IX)*UU(2,IX+1)  
     .                                + FF(1,IX)               
               UU(2,IX) = EE(2,1,IX)*UU(1,IX+1) + EE(2,2,IX)*UU(2,IX+1)  
     .                                    + FF(2,IX)               
               VZN(IX,J) = AAA(IX)*VZN(IX+1,J) + BBB(IX)                 
        !#! Marco
!         if (m .eq. 0 .and. n .eq. 0 .and. ix .eq. 50) then
!          write(*,*) 'semimpl', EE(1,1,IX),
!     .                EE(1,2,ix), EE(2,1,ix), EE(2,2,ix),
!     .                FF(1,ix),FF(2,ix),
!     .                AAA(ix), BBB(ix)
!         endif
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
                                                                                
810     DO 351 IX=0,LX                                                   
               VRN(IX,J1) = CONJG(VRN(IX,J))                                    
351     CONTINUE                                                          
               DO 352 IX=-1,LX                                                  
               VTN(IX,J1) = CONJG(VTN(IX,J))                                    
               VZN(IX,J1) = CONJG(VZN(IX,J))                                    
352     CONTINUE                                                          
15       CONTINUE                                                          
                                                                        
!#! Marco
!        write(*,*),'intermedio0 ,it= ',it
!     .             ,'btn0=',btn(50,j0),' bt0=',bt(50,j0)
!     .             ,'brn0=',brn(50,j0),' br0=',br(50,j0)
!     .             ,'vtn0=',vtn(50,j0),' vt0=',vt(50,j0)
!     .             ,'brn=',brn(100,74),' br=',br(100,74)
!     .             ,'bzn=',bzn(100,74),' bz=',bz(100,74)

!          write(*,*)  'medio2= ',bz(90,j0)

C ............................  B  .....................................
!#! Marco, seconda legge di Faraday
       DO 28 IM=0,MY                                                    
       M = MM(IM)                                                       
       DO 28 N=NA1(M),NE(M)                                             
          SUM0=(0.0,0.0)
          SUM1=(0.0,0.0)
          SUM2=(0.0,0.0)
          SUM3=(0.0,0.0)
          SUM8=(0.0,0.0)
       jmn = janz(m,n)
!       if (real(vtn(50,jmn)) .gt. 10.d0) then
!        write(*,*) 'fotte m,n ',m , n, vtn(50,jmn),
!     .            vt(50,jmn), btn(50,jmn), bt(50,jmn)
!       endif
       if (m.eq.0) jmn = -jmn ! serve altrimenti ho jmn negativi!!!
!#! Marco, voglio anche m,n come variabili privare
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(IAM,js,IX,J,J1,m,n)
!$OMP& SCHEDULE(STATIC)
!$OMP& REDUCTION(+:SUM0,SUM1,SUM2,SUM3,SUM8)
       DO 30 js=0,convnj(jmn)-1
       IAM=OMP_GET_THREAD_NUM()
       j = convj(jmn,js)
       j1 = convj1(jmn,js)
       DO 32 IX=0,LX-1                                                  
       SUM3(IX) = SUM3(IX)+(VT(IX,J)+VTN(IX,J))*BZN(IX,J1)              
     . -(VZ(IX,J)+VZN(IX,J))*BTN(IX,J1)                                 
   32     CONTINUE                                                          
       DO 33 IX=0,LX                                                    
       SUM0(IX) = SUM0(IX)+BRN(IX,J)*                                   
     . (VT(IX-1,J1)+VTN(IX-1,J1)+VT(IX,J1)+VTN(IX,J1))                  
       SUM1(IX)=SUM1(IX)+(VR(IX,J)+VRN(IX,J))*(BTN(IX-1,J1)+BTN(IX,J1)) 
       SUM2(IX) = SUM2(IX)+BRN(IX,J)*                                   
     . (VZ(IX-1,J1)+VZN(IX-1,J1)+VZ(IX,J1)+VZN(IX,J1))
       SUM8(IX) = SUM8(IX)+(VR(IX,J)+VRN(IX,J))*                        
     . (BZN(IX-1,J1)+BZN(IX,J1))                                        
!#! Marco
!       if (m .eq. 0 .and. n .eq. 0 .and. ix .eq. 50
!     .       .and. j .eq. 0 .and. j1 .eq. 0) then
!!!        write(*,*)  'an= ',js,j,j1,
!        write(*,*)  'an2= ',js,j,j1
!     .            ,sum0(ix), sum1(ix),brn(ix,j)
!     .            ,vt(ix-1,j1),vtn(ix-1,j1),vt(ix,j1),vtn(ix,j1)
!!!     .           br(ix,j), vr(ix,j),
!!!     .           vt(ix,j1), bt(ix,j1)
!       endif
!!#! Marco
!       if (ix .eq. 50 .and. js .eq. 0) then
!       write(*,*)  iam,'ancora piu Far2= ',js,j,j1,
!     .           brn(ix,j),
!     .           vrn(ix,j),
!     .           vtn(ix,j), vzn(ix,j)
!       endif
33     CONTINUE                                                          
30     CONTINUE                                                          
!$OMP END PARALLEL DO

!          write(*,*)  'dopo seconda faraday,it= ',it,
!!     .          'brn=',brn(50,j0),' br=',br(50,j0),
!!     .          'btn=',btn(50,j0),' bt=',bt(50,j0),
!     .          'bzn=',bzn(90,j0),' bz=',bz(90,j0)

C ACHTUNG MANCHE SUMMEN (ABL.) AUCH AM RAND NOTWENDIG                   
       DO 34 IX=0,LX-1                                                  
       DSUM2(IX)=(SUM1(IX+1)-SUM1(IX)-SUM0(IX+1)+SUM0(IX))*ZS22(IX)*0.5 
       DSUM3(IX) =((SUM2(IX+1)-SUM8(IX+1))*X(IX+1)-(SUM2(IX)-SUM8(IX))  
     . *X(IX))*ZS22(IX)/X2(IX)*0.5                                      
!#! Marco
!       if (m .eq. 0 .and. n .eq. 0 .and. ix .eq. 50) then
!        write(*,*) 'guardo dsum2,seconda faraday ',
!     .           sum1(ix),sum0(ix),dsum2(ix)
!       endif
34     CONTINUE                                                          
       J = JANZ(MM(IM),N)                                               
       J1 = JANZ(-MM(IM),-N)                                            
       CFAK1 = N/RR*I                                                   
       DO 35 IX=0,LX-1                                                  
       CFAK = M/X2(IX)*I                                                
       BT(IX,J) = BT(IX,J)+DT2*(+CFAK1*SUM3(IX)-DSUM2(IX))              
       BZ(IX,J) = BZ(IX,J)+DT2*(+DSUM3(IX)-CFAK*SUM3(IX))               
35     CONTINUE                                                          
!       if (j eq j0) then
!        write(*,*),'dopo dt2 ,it= ',it,
!     .             'btn=',btn(50,j0),' bt=',bt(50,j0)
!       endif

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
!#! Marco
!!       if (j eq j0) then
!        write(*,*),'intermedio1 ,it= ',it,
!     .             'btn=',btn(50,j0),' bt=',bt(50,j0)
!     .             ,'vtn=',vtn(50,j0),' vt=',vt(50,j0)
!!       endif
!          write(*,*)  'medio3= ',bz(90,j0)

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
!          write(*,*)  'medio4= ',bz(90,j0)
C  BR AUS BT +BZ                                                        
74     BR(0,J) =(0.0,0.0)                                               
       DO 150 IX=1,LX-1                                                 
       BR(IX,J)=(I/ZS22(IX-1)*(-BT(IX-1,J)*M-BZ(IX-1,J)*N*X2(IX-1)      
     . /RR)+BR(IX-1,J)*X(IX-1))/X(IX)                                   
150     CONTINUE                                                          
       IF(M.NE.1) GO TO 73                                               
C      BR(0,J) = 0.5*(BR(1,J)-BT(1,J))                                  
       BR(0,J) = (4.0*BR(1,J)-BR(2,J))*0.333333333                      
       BT(-1,J) = -BT(0,J)+2.0*I*BR(0,J)                                
73     IF(J.EQ.0) GO TO 82                                               
       DO 40 IX=0,LX                                                    
       BR(IX,J1) = CONJG(BR(IX,J))                                      
40     CONTINUE                                                          
       DO 402 IX=-1,LX                                                  
       BT(IX,J1) = CONJG(BT(IX,J))                                      
       BZ(IX,J1) = CONJG(BZ(IX,J))                                      
402     CONTINUE                                                          
82     CONTINUE                                                          
28     CONTINUE                                                          
!#! Marco
!        write(*,*),'intermedio2 ,it= ',it
!     .             'btn=',btn(50,j0),' bt=',bt(50,j0)
!     .             ,'brn=',brn(100,74),' br=',br(100,74)
!     .             ,'bzn=',bzn(100,74),' bz=',bz(100,74)

!          write(*,*)  'fine= ',bz(90,j0)
!       WRITE(6,*) 'dopo28   0.5*(REAL(BT(LX-1,0))+REAL(BT(LX,0))) = ',
!     .        0.5*(REAL(BT(LX-1,0))+REAL(BT(LX,0)))
                                                                        
       DO 41 IM=0,MY                                                    
       M = MM(IM)                                                       
       DO 41 N=NA1(M),NE(M)                                             
       J = JANZ(M,N)                                                    
       J1= JANZ(-M,-N)                                                  
         DO 43 IX=0,LX-1                                                
              R3(IX) = BZ(IX,J) - DT*ZDR*(ETAA(IX+1)-ETAA(IX))*0.5*       
     .                 I*N/RR*(BR(IX,J) + BR(IX+1,J))                   
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
     .             - A13*ZDR*X(1)/X2(0)                                 
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
     .                     - A13/X2(IX)                                 
                                                                        
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
     .                           (BETJ(IX) + GAMJ(IX)*AAA(IX-1))        
902     CONTINUE                                                       
                                                                        
!         RRAD = X2(LX-1)/(2. - X2(LX-1))
! Daniele, 25/11/2008
! caso generale per VRLX =/ 0
         RRAD = (2.*ETAA(LX)*X2(LX-1) + VRLX/LX)/
     .          (2.*ETAA(LX)*X2(LX) - VRLX/LX)
! Daniele, 23/06/2009
! questo è per imporre il campo magnetico radiale al bordo=0.
         UU(1,LX) = (0.0,0.0)

!         IF (J.EQ.79) UU(1,LX) = (0.0,1.5E-2)

         IF (J.EQ.0) THEN                                               
! Daniele, 18/08/2009
! caso generale per VRLX =/ 0 e BR =/0
              UU(2,LX-1) = (EE(2,1,LX-1)*UU(1,LX)+
     .               (EE(2,2,LX-1)*E0*(2./LX))
     .             /(2.*ETAA(LX)*X2(LX) - VRLX/LX) + FF(2,LX-1))
     .             /(1. - EE(2,2,LX-1)*RRAD)          
              UU(2,LX) = UU(2,LX-1)*RRAD +
     .               E0*(2./LX)/(2.*ETAA(LX)*X2(LX) - VRLX/LX)
!              write(*,*) 'in the beg',E0,2./LX,UU(1,LX)

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

C      .................                                                
         ELSE                                                           

!#! Marco, applico MP irrotazionali
         call apply_mp

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
     .                                         + FF(1,IX)               
              UU(2,IX) = EE(2,1,IX)*UU(1,IX+1) + EE(2,2,IX)*UU(2,IX+1)  
     .                                         + FF(2,IX)               
              BZN(IX,J) = AAA(IX)*BZN(IX+1,J) + BBB(IX)                 
907     CONTINUE                                                       
         DO 908 IX=0,LX                                                 
              BRN(IX,J) = UU(1,IX)                                      
              BTN(IX,J) = UU(2,IX)                                      
908     CONTINUE                                                       

!         if (j .eq. 0) then 
!          write(*,*),'pre-solenoidalita2 ,it= ',it,m,n,j
!     .             ,'btn0=',btn(50,j0),' bt0=',bt(50,j0)
!     .             ,'brn0=',brn(50,j0),' br0=',br(50,j0)
!!         endif
!     .             ,'vtn=',vtn(50,j0),' vt=',vt(50,j0)
!     .             ,'brn=',brn(100,74),' br=',br(100,74)
!     .             ,'bzn=',bzn(100,74),' bz=',bz(100,74)
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
909      IF (M.NE.1) GO TO 910                                          
         BTN(-1,J) = - BTN(0,J) + 2.*I*BRN(0,J)                         
910     CONTINUE                                                       
                                                                        
C  BR AUS BT +BZ                                                        
       BRN(0,J) =(0.0,0.0)                                              
       
!#! Marco, calcolo anche il punto LX
       DO 151 IX=1,LX-1                                                 
!       DO 151 IX=1,LX
       BRN(IX,J)=(I/ZS22(IX-1)*(-BTN(IX-1,J)*M-BZN(IX-1,J)*N*X2(IX-1)   
     . /RR)+BRN(IX-1,J)*X(IX-1))/X(IX)                                  
151     CONTINUE                                                          
      IF(M.NE.1) GO TO 75                                               
       BRN(0,J) = (4.0*BRN(1,J)-BRN(2,J))*0.33333333333                 
       BTN(-1,J) = -BTN(0,J)+2.0*I*BRN(0,J)                             
                                                                        
75     IF(J.EQ.0) GO TO 83                                               
       DO 45 IX=0,LX                                                    
       BRN(IX,J1) = CONJG(BRN(IX,J))                                    
45     CONTINUE                                                          
       DO 404 IX=-1,LX                                                  
       BTN(IX,J1) = CONJG(BTN(IX,J))                                    
       BZN(IX,J1) = CONJG(BZN(IX,J))                                    
404   CONTINUE                                                          
83    CONTINUE                                                          
41     CONTINUE                                                          
!#! Marco
!        write(*,*),'solenoidalita2 ,it= ',it
!     .             ,'btn0=',btn(50,j0),' bt0=',bt(50,j0)
!     .             ,'brn0=',brn(50,j0),' br0=',br(50,j0)
!     .             ,'vtn0=',vtn(50,j0),' vt0=',vt(50,j0)
!     .             ,'brn=',brn(100,74),' br=',br(100,74)
!     .             ,'bzn=',bzn(100,74),' bz=',bz(100,74)


!#! Marco, azzero d'autorità (per ora), la parte immaginaria delle
!componenti poloidale e toroidale del campo magnetico di equilibrio.
        do ix=-1,lx 
           bt(ix,j0) = dcmplx(real(bt(ix,j0)), 0.d0)
           bz(ix,j0) = dcmplx(real(bz(ix,j0)), 0.d0)
           btn(ix,j0) = dcmplx(real(btn(ix,j0)), 0.d0)
           bzn(ix,j0) = dcmplx(real(bzn(ix,j0)), 0.d0)
!#! Marco, 10/10/2014 azzero anche la parte immaginaria delle velocità
!(che nelle simulazioni è dell'ordine di  10^-19
           vrn(ix,j0) = dcmplx(real(vrn(ix,j0)), 0.d0)
           vtn(ix,j0) = dcmplx(real(vtn(ix,j0)), 0.d0)
           vzn(ix,j0) = dcmplx(real(vzn(ix,j0)), 0.d0)
        enddo

!#! Marco
!        write(*,*)  'fine2,it= ',it,
!     .             'brn=',brn(50,j0),' br=',br(50,j0),
!     .             'btn=',btn(50,j0),' bt=',bt(50,j0),
!     .             'bzn=',bzn(50,j0),' bz=',bz(50,j0)
                                                                        
C..................  J-BESTIMMUNG ..................................... 
                                                                        
       DO 425 IM=0,MY                                                   
       M = MM(IM)                                                       
       DO 425 N=NA1(M),NE(M)                                            
       FAK1 = N/RR                                                      
       J = JANZ(M,N)                                                    
       J1 = JANZ(-M,-N)                                                 
       DO 426 IX=0,LX-1                                                 
       JR(IX,J) = I*(M/X2(IX)*BZN(IX,J)-FAK1*BTN(IX,J))                 
426     CONTINUE                                                          
!       DO 427 IX=1,LX-1                                                 
       DO 427 IX=1,LX
       JT(IX,J) = I*FAK1*BRN(IX,J)-2.0*ZS2(IX)*(BZN(IX,J)-BZN(IX-1,J))  
       JZ(IX,J) = 2.0*ZS2(IX)/X(IX)*(X2(IX)*BTN(IX,J)-X2(IX-1)          
     . *BTN(IX-1,J))-I*M/X(IX)*BRN(IX,J)                                
427    CONTINUE                                                          
!       JT(LX,J) = (0.0,0.0)                                             
!       JZ(LX,J) = (0.0,0.0)                                             
       JT(0,J) =(0.0,0.0)                                               
       JZ(0,J) =(0.0,0.0)                                               
       IF(M.NE.0) GO TO 422                                             
!       JT(LX,J) = -2.*ZS2(LX)*(BZN(LX,J) - BZN(LX-1,J))                 
!       JZ(LX,J) =  2.*ZS2(LX)*((2. - X2(LX-1))*BTN(LX,J) -              
!     &                                          X2(LX-1)*BTN(LX-1,J))   
        JZ(0,J) = 2.*BTN(0,J)/X2(0)                                     
422    IF(M.NE.1) GO TO 440                                              
       JT(0,J) = JT(1,J)                                                
440    IF(J.EQ.0) GO TO 425                                              
       DO 439 IX=0,LX                                                   
       JR(IX,J1) = CONJG(JR(IX,J))                                      
439    CONTINUE                                                          
       DO 419 IX=0,LX                                                   
       JT(IX,J1) = CONJG(JT(IX,J))                                      
       JZ(IX,J1) = CONJG(JZ(IX,J))                                      
419    CONTINUE                                                          
425    CONTINUE                                                          
C                                                                       
                                                                        
! ............................DENSITY ..................................
!                                                                       
!     DO 720 IM=0,MY                                                    
!      M = MM(IM)                                                       
!      DO 720 N=NA1(M),NE(M)                                            
!                                                                       
!       DO 721 IX = 0,LX                                                
!       SUM0(IX) = (0.0,0.0)                                            
!       SUM1(IX) = (0.0,0.0)                                            
!721     CONTINUE                                                       
!                                                                       
!       DO 722 MS = -MM(MY)+M,MM(MY)                                    
!      IF(MWERT(MS).EQ.0.OR.MWERT(M-MS).EQ.0) GO TO 722                 
!        DO 723 NS = NA(MS),NE(MS)                                      
!      IF(N-NS.LT.NA(M-MS).OR.N-NS.GT.NE(M-MS)) GO TO 723               
!        J = JANZ(MS,NS)                                                
!        J1= JANZ(M-MS,N-NS)                                            
!          DO 724 IX=1,LX-1                                             
!                                                                       
!          VTIXJ= ((VT(IX-1,J)+VTN(IX-1,J))/2. +                        
!    .             (VT(IX,J)  +VTN(IX,J)  )/2.  )/2.                    
!          VZIXJ= ((VZ(IX-1,J)+VZN(IX-1,J))/2. +                        
!    .             (VZ(IX,J)  +VZN(IX,J)  )/2.  )/2.                    
!                                                                       
!     DENSITY :                                                         
!          S0 = - ZDR * 0.5/ X(IX)                                      
!          S1 = X(IX+1)*DENN(IX+1,J1)*(VR(IX+1,J)+VRN(IX+1,J))/2.       
!          S2 = X(IX-1)*DENN(IX-1,J1)*(VR(IX-1,J)+VRN(IX-1,J))/2.       
!          S3 = DENN(IX,J1)                                             
!          S4 = I * M * VTIXJ  / X(IX)                                  
!          S5 = I * N * VZIXJ   / RR                                    
!          SUM0(IX) = S0 *(S1-S2) - S3 *(S4+S5) + SUM0(IX)              
!                                                                       
!    PRESSURE :                                                         
!          S0 = - ZDR * 0.5 / X(IX)                                     
!          S1 = X(IX+1)*PRNN(IX+1,J1)*(VR(IX+1,J)+VRN(IX+1,J))/2.       
!          S2 = X(IX-1)*PRNN(IX-1,J1)*(VR(IX-1,J)+VRN(IX-1,J))/2.       
!          S3 = PRNN(IX,J1)                                             
!          S4 = I * M * VTIXJ  / X(IX)                                  
!          S5 = I * N * VZIXJ   / RR                                    
!          S6 = -2./3. * PRNN(IX,J1) * (                                
!    .          1./X(IX) * (                                            
!    .     X(IX+1)* (VR(IX+1,J)+VRN(IX+1,J))/2.  -                      
!    .     X(IX-1)* (VR(IX-1,J)+VRN(IX-1,J))/2.    ) *0.5 *ZDR          
!    .       + I * MS /X(IX) * VTIXJ + I * NS /RR * VZIXJ    )          
!          SUM1(IX) = S0 *(S1-S2) - S3 *(S4+S5) + S6 + SUM1(IX)         
!                                                                       
!724       CONTINUE                                                     
!723       CONTINUE                                                     
!722       CONTINUE                                                     
!                                                                       
!          J = JANZ(MM(IM),N)                                           
!          J1= JANZ(-MM(IM),-N)                                         
!          DO 725 IX=1,LX-1                                             
!                                                                       
!          AMX = 0.2                                                    
!          DEN(IX,J)=AMX*(DE(IX+1,J)+DE(IX-1,J))/2.+ (1.-AMX)*DE(IX,J)  
!    .                +  DT * SUM0(IX)                                  
!                                                                       
!          PRN(IX,J)=AMX*(PR(IX+1,J)+PR(IX-1,J))/2.+ (1.-AMX)*PR(IX,J)  
!    .                +  DT * SUM1(IX)                                  
!                                                                       
!          DEN(IX,J1)= CONJG(DEN(IX,J))                                 
!          PRN(IX,J1)= CONJG(PRN(IX,J))                                 
!725       CONTINUE                                                     
!          IF ( MM(IM) .EQ. 0 ) THEN                                    
!          DEN(0,J)  = DEN(1,J)                                         
!          DEN(0,J1) = CONJG(DEN(0,J))                                  
!          PRN(0,J)  = PRN (1,J)                                        
!          PRN(0,J1) = CONJG(PRN(0,J))                                  
!          ELSE                                                         
!          DEN(0,J)  = 0.0                                              
!          DEN(0,J1) = CONJG(DEN(0,J))                                  
!          PRN(0,J)  = 0.0                                              
!          PRN(0,J1) = CONJG(PRN(0,J))                                  
!          ENDIF                                                        
!          IF( J .EQ. 0 )THEN                                           
!          DEN(LX,J) = DEN(LX-1,J)                                      
!          ELSE                                                         
!          DEN(LX,J) = 0.0                                              
!          ENDIF                                                        
!          DEN(LX,J1)=CONJG(DEN(LX,J))                                  
!         PRN(LX,J) = PRN(LX-1,J)                                       
!         PRN(LX,J) = PR(LX,J)                                          
!         PRN(LX,J1)=CONJG(PRN(LX,J))                                   
!720       CONTINUE                                                     
!                                                                       
!......................... correzione del flusso toroidale:
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

        ENDIF

! daniele, 11/05/2007
        IF( TORFLX .GT. TORF1*(1.d0+0.005d0) ) THEN

        TORF = 0.0
        DO 989 IX=0,LX-1
        TORF=TORF+ 2.* X2(IX)*REAL(BZ(IX,J0))/ZDR
989     CONTINUE
        FTOR = TORF1 / TORF
        WRITE(6,*) '  guadagnato 0.5% del flusso rispetto a it= 0 '
        WRITE(6,*) '  correzione, it = ',it
        call flush(6)

        ENDIF
       ENDIF

!!#! Marco, 07/01/2015, rotating_mp, intanto le metto qui, prima della
!reinizializzazione per il timestep successivo
         if (irrot .eqv. .true. .and. mp_perc(0) .ne. 0.d0
     .       .and. br_pert_freq /= 0.d0) then
!#! REMARK: tolgo e aggiungo a B?N

!#! inizioparte2
           do idx=0, size(n_mp) -1

            call apply_irrot_bc(mp_perc(idx),amp_a(idx),m_mp(idx)
     .                         ,n_mp(idx),phmp(idx),bbr,bbt,bbz)
!#! Marco, subtract irrotational perturbation 
            do ix=0, lx 
             brn(ix,j_mp(idx)) = brn(ix,j_mp(idx)) - bbr(ix)
             brn(ix,-j_mp(idx)) = conjg(brn(ix,j_mp(idx)))
            enddo
            do ix=-1, lx 
             btn(ix,j_mp(idx)) = btn(ix,j_mp(idx)) - bbt(ix)
             bzn(ix,j_mp(idx)) = bzn(ix,j_mp(idx)) - bbz(ix)
             btn(ix,-j_mp(idx)) = conjg(btn(ix,j_mp(idx)))
             bzn(ix,-j_mp(idx)) = conjg(bzn(ix,j_mp(idx)))
            enddo

!            if (j_mp(idx) .eq. 73) then
!             write(*,*) 'ph0',it,phmp(idx),amp_a(idx)
!             write(*,*) 'br0',brn(lx,73),abs(brn(lx,73))
!             write(*,*) 'add0',bbr(lx),abs(brn(lx,73))
!            endif

!#! Marco, update phase
            phmp(idx) = phmp(idx) + br_pert_freq*dt
            if (phmp(idx) .gt. 2.d0*pi) then
             phmp(idx) = phmp(idx) - 2.d0*pi
            endif
           
!#! Marco, add irrotational perturbation 
            call apply_irrot_bc(mp_perc(idx),amp_a(idx),m_mp(idx)
     .                          ,n_mp(idx),phmp(idx),bbr,bbt,bbz)
            do ix=0, lx 
             brn(ix,j_mp(idx)) = brn(ix,j_mp(idx)) + bbr(ix)
             brn(ix,-j_mp(idx)) = conjg(brn(ix,j_mp(idx)))
            enddo
            do ix=-1, lx 
             btn(ix,j_mp(idx)) = btn(ix,j_mp(idx)) + bbt(ix)
             bzn(ix,j_mp(idx)) = bzn(ix,j_mp(idx)) + bbz(ix)
             btn(ix,-j_mp(idx)) = conjg(btn(ix,j_mp(idx)))
             bzn(ix,-j_mp(idx)) = conjg(bzn(ix,j_mp(idx)))
            enddo
!            if (j_mp(idx) .eq. 73) then
!             write(*,*) 'ph1',it,phmp(idx),amp_a(idx)
!             write(*,*) 'br1',brn(lx,73),abs(brn(lx,73))
!             write(*,*) 'add1',bbr(lx),abs(brn(lx,73))
!            endif
     
!#! Marco. Noto che il modulo di br(lx,jmp) varia nel tempo. Invento un
!fattore di correzione.
            corr_fact = mp_perc(idx) / abs(brn(lx,j_mp(idx)))
            write(*,*) 'DIAG_corr_fact=', corr_fact
            brn(:,j_mp(idx)) = brn(:,j_mp(idx)) * corr_fact
            btn(:,j_mp(idx)) = btn(:,j_mp(idx)) * corr_fact
            bzn(:,j_mp(idx)) = bzn(:,j_mp(idx)) * corr_fact
!          write(*,*) bbr
!          write(*,*) '*****'
!          write(*,*) bbz
!          write(*,*) '*****'

!#! fineparte2
           enddo
         endif

!       write(*,*) 'alla fine',brn(50,0),br(50,0)
!     .            ,btn(50,0),bt(50,0)
!     .            ,vtn(50,0),vt(50,0)

C......................... REINIZIALIZZAZIONE :                         
992    DO 77 IM=0,MY                                                    
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
403     CONTINUE 
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
405       CONTINUE                                                         
77       CONTINUE                                                         

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
!#! Marco, 06/08/2014 I write the whole complex number
      WRITE(6,*) '   BZ(LX-1,J0) =',BZ(LX-1,J0),
     .        '   (BZ(LX,J0) = ', BZ(LX,J0)
       WRITE(6,*) '   BT(LX-1,J0) =',BT(LX-1,J0),
     .        '   BT(LX,J0) = ', BT(LX,J0)
       WRITE(6,*) '   Bt(LX,73) =',BT(LX,73),
     .        '   Bz(LX,73) = ', Bz(LX,73),
     .            'Br(LX,73) = ', BR(lx,73)
!       WRITE(6,*) '   0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0))) = ',
!     .        0.5*(REAL(BT(LX-1,J0))+REAL(BT(LX,J0)))
      WRITE(6,*) '    flusso toroidale= ', TORFLX,'  ftor=',FTOR
      WRITE(6,*) '   TH =',YTH, '   F = ', YF(ITP)
      WRITE(6,*) '    ITP=', ITP
      WRITE(6,*) '    time=', it*dt
      if (br_pert_freq /= 0.d0) then
       write(*,*) 'phmp = ',phmp
      endif
      call flush(6)
                                                                        
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
                                                
         IF( IT .EQ. ITEND ) then
!          WRITE(IBAUS) ZEIT
!          write(ibaus) int(lx),int(iz)
          WRITE(IBAUS) ZEIT,lx,iz
!          write(*,*) 'scrivo in end',zeit,lx,iz
         endif
                                                
       WRITE(IBOUT1)  ZEIT,lx,iz
!       write(ibout1)  int(lx),int(iz)
!          write(*,*) 'scrivo in all',zeit,lx,iz
!#! Marco 06/08/2014
!#! I write the whole complex number, not only the real or imaginary
!part in the fixed-phases assumptions
       DO 252 J=0,IZ                            
          IF ( IT .EQ. ITEND ) then
           WRITE(IBAUS)   VT(-1,j), VZ(-1,j), BT(-1,j), BZ(-1,j) 
          endif
          WRITE(IBOUT1)  VT(-1,j), VZ(-1,j), BT(-1,j), BZ(-1,j)
                                                
       DO 252 IX=0,LX                           
                                                                        
          IF (IT .EQ. ITEND) then
           WRITE(IBAUS) VR(ix,j), VT(ix,j), VZ(ix,j),
     .                  BR(ix,j), BT(ix,j), BZ(ix,j)
          endif
          WRITE(IBOUT1) VR(ix,j), VT(ix,j), VZ(ix,j), 
     .                  BR(ix,j), BT(ix,j), BZ(ix,j) 
                                                                        
 252   CONTINUE                                                         

       call flush(IBAUS)
       call flush(IBOUT1)
 

 117   CONTINUE                                                         
   50 CONTINUE                                                          

       close(unit=ibaus)
      close(unit=ibout1)
5000  continue

C................ MAIN LOOP ENDS HERE ................................  
1111         write(6,*) ' zeit= ',zeit,'  stato finale:'
         WRITE(6,*) 'x= ',X2(LX-4),' bz ',BZ(LX-4,0),' bt ',BT(LX-4,0)
         WRITE(6,*) 'x= ',X2(LX-3),' bz ',BZ(LX-3,0),' bt ',BT(LX-3,0)
         WRITE(6,*) 'x= ',X2(LX-2),' bz ',BZ(LX-2,0),' bt ',BT(LX-2,0)
         WRITE(6,*) 'x= ',X2(50),' bz ',BZ(50,0),' bt ',BT(50,0)
         WRITE(6,*) 'x= ',X2(50),' br ',BR(LX,73),' bt ',BT(LX,73)
        
         call flush(6)

       call dealloc_arrays

  100 FORMAT(6E12.4)                                                    
  101 FORMAT(I4)                                                        
  102 FORMAT(10E12.4)                                                   
  103 FORMAT(E12.4,2I4)                                                 
  550 FORMAT(1X,' ITEND,ITP=',2I10,' ZEIT=',1PE12.4,' RS=',E12.4)       
       END                                                              
