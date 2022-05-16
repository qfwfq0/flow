
        program pppoincare_par
        
        use modpoinc

        use defs
	

	character gran*20,t1*30,t2*20,t3*20,t4*20
	character titp*20,tlth*10,tlphi*10,cmd*100

!	integer, parameter :: lth=32
!	integer, parameter :: lphi=64
	integer(wp) :: lth,lphi
	integer(4) :: ly2,lth2,lphi2
	integer(wp) :: iphimod
!!$OMP THREADPRIVATE(lth,lphi)
!      OMP variables
	integer,save :: thr_tot,thr_num
!$OMP THREADPRIVATE(thr_tot,thr_num)
	integer :: omp_get_num_threads,omp_get_thread_num
!	integer,save :: im,n,m
!$OMP THREADPRIVATE(im,n,m,jmn)
	integer,save :: iphi,ith!,ix
!$OMP THREADPRIVATE(iphi,ith,ix)
        real(wp),save :: fmns,fmnc,fons,fonc,foo,ieli,rnd
!$OMP THREADPRIVATE(fmns,fmnc,fons,fonc,foo,ieli,rnd)
        real(wp),save :: aphi,ath
        real(wp) :: ppi
!$OMP THREADPRIVATE(aphi,ath)
        real(wp),save :: abt,abz,fact,dom_fact
        real(wp), save :: abrix
!$OMP THREADPRIVATE(abrix,abt,abz,fact,dom_fact)
        real(wp),save :: nmode!, r0
!$OMP THREADPRIVATE(nmode, r0)
        real(wp), allocatable, dimension(:) :: vecr
        real(wp), allocatable, dimension(:)  :: vecth
        real(wp), allocatable, dimension(:)  :: vecphi
!        complex(8),allocatable,dimension(:,:) :: bbr
!        complex(8),allocatable,dimension(:,:) :: bbth
!        complex(8),allocatable,dimension(:,:) :: bbphi
        real(wp),allocatable,dimension(:,:,:) :: matr
        real(wp),allocatable,dimension(:,:,:) :: matth
        real(wp),allocatable,dimension(:,:,:) :: matphi
        real(wp),allocatable,dimension(:,:,:) :: matb1
        real(wp),allocatable,dimension(:,:,:) :: matb2
        real(wp),allocatable,dimension(:,:,:) :: matb3
!        integer :: uft=41
	integer(4) :: bcond(6)
        integer(4) :: izzz,lyy
        real(8) :: zeitt
        integer(wp) :: section
        logical :: result
	logical :: tor=.false.
        logical,save :: mode_dump
        logical,save :: shsh
        logical,save :: dominant_dump=.false.
        logical,save :: random_phases=.false.
!$OMP THREADPRIVATE(mode_dump,dominant_dump,random_phases)
	integer, dimension (9):: nfiles
        integer :: inlun, nmodes,  d3,len_itp
        real(wp) :: inteli
        integer :: modenm(0:1)
        integer, dimension(0:10) :: iitp
        character(300) :: path,inputfile,tfile
        character(30) :: ditp, heli, dnth, dnphi
        character(len=255) :: nmtfile

        namelist /nmt/ iitp,inteli,lth,lphi,shsh,mode_dump
!#! Marco, load default values of the spc2 namelist
        iitp(0:10) = 0
        inteli = 0
        lth = 32
        lphi = 64
        mode_dump = .false.
        shsh = .false.

        namelist /spc2/ lx,my,iz,i2d,n2d,m2d
     $   ,fasi_sganciate,irrot

!#! I define a new namelist containing relevant infos about spectra,
!equilibrium, etc
      namelist /spc3/ qmp,br_pert_freq,mp_perc,mp_cmplx
     $   ,ph_mp,m_mp,n_mp
     $   ,mm,nanf,nz,nzcon,m0_momsour
     $   ,m_mpold,n_mpold,ampold,pitol,tok_rot,tkp,ipr
     $   ,chi,beta0,theta0,alpha,q01,aa1,bb1,rr,dt
     $   ,itend,ib,itor,eta0,alet,beet,rho0,gaet,mue,almu,bemu
     $   ,lambda,iband,iband1,inopert
     $   ,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
     $   ,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20
     $   ,tau_growth,tau_decrease,time_passage
     $   ,mpexp,mptanh,time_medium

!#! preliminary I-O operations
        call getcwd(prefix)
!        unit_all = trim(trim(prefix)//"/dat/"
!     .           //"specyl_"//trim(string)//"_all.dat")
  
  !#! Marco, load default values of the spc2 namelist
        call default_spc2 
   
        spc2file = trim(trim(prefix)//"/for/1/spc2.in")
        spc3file = trim(trim(prefix)//"/for/1/spc3.in")
  
        write(*,*) 'spc2in',spc2file
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
	write(*,*) 'check_alloc',my,iz,lx
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
  !      write(*,*) 'open3',ios
        read(nunitd,spc3,iostat=ios)
  !      write(*,*) 'read3',ios
        if (ios /= 0) then
         write(*,*) 'could not read namelist spc3.in'
         write(*,*) spc3file
         stop
  
        else
         close(nunitd)
        endif

        nmtfile = trim(trim(prefix)//"/for/nmt.in")

        open(unit=uft, file=trim(nmtfile)
     .     ,status='unknown', iostat=ios)
        if (ios /= 0) then
         write(*,*) 'could not open namelist nmt.in'
         stop
        endif
        read(uft,nmt,iostat=ios)
        if (ios /= 0) then
         write(*,*) 'could not read namelist nmt.in'
        else
         close(nunitd)
        endif

        bcond(1) = 7
        bcond(2) = 0
        bcond(3) = 1
        bcond(4) = 1
        bcond(5) = 1
        bcond(6) = 1

        d3 = 1
        r0 = 4.d0
        write(*,*) 'SH=',shsh
        write(*,*) 'inteli=',inteli
!        read( heli, '(i10)' ) inteli
!        read( dnth, '(i10)' ) lth
!        read( dnphi, '(i10)' ) lphi
        
!#! Marco, opening and reading the file containing the magnetic field
!profiles
        section = 1
        it = 0

        do q=0,10 
!#! Marco, apro il primo file
        unit_all = trim(trim(prefix)//"/dat/"
     .        //"specyl_"//trim(aint2char(section))//"_all.dat")
        open(unit=ibein,file=trim(unit_all),action='read'
     .          ,form='unformatted')

        if (iitp(q) .gt. 0) then
         do while (it .le. iitp(q))

           read(ibein,iostat=ierr) zeitt,lyy,izzz

!           if (abs(mod(zeitt,10000.)) .lt. 1.d0) then
!            write(*,'(a,f14.7,a,i4,a,i4)') 
!     .         'zeitt=',zeitt,' lyy=',lyy,'  izzz=',izzz
!           endif
            write(*,'(a,f14.7,a,i0,a,i0,a,i0,a,i0)') 
     .         'zeitt=',zeitt,' lyy=',lyy,' izzz=',izzz,'  it=',it,'
     .           err=',ierr
           if (ierr .lt. 0) then !end of file reached
            close (ibein)
            it = it - 1
            section = section + 1
            write(*,*) 'ho aperto lo specyl_?.dat',section
            unit_all = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(aint2char(section))//"_all.dat")
            open(unit=ibein,file=trim(unit_all),action='read'
     .        ,form='unformatted')
           read(ibein,iostat=ierr) zeitt,lyy,izzz
!            if (mod(zeitt,100.) .eq. 0.d0) then
!             write(*,*) 'file successivo, section=',section
!             write(*,'(a,f9.5,a,i4,a,i4)') 
!     .         'zeitt=',zeitt,' lyy=',lyy,'  izzz=',iz
!            endif
           endif

         call flush(6)
         if( izzz .ne. iz ) then
          write(6,*) 'numero diverso di modi elaborati',iz,izzz
          stop
         endif
         call flush(6)
       
         write(*,*) 'fin qui 0'
         
         do j=0,izzz
          read(ibein) vtold(-1,j), vzold(-1,j), btold(-1,j), bzold(-1,j)
           do ix=0,lyy                                                   
            read(ibein) vrold(ix,j), vtold(ix,j), vzold(ix,j),
     .                 brold(ix,j), btold(ix,j), bzold(ix,j) 
           enddo
         enddo  

         if (it .eq. iitp(q)) then 
          write(*,*) 'sto leggendo, ma ho trovato quindi...'
          write(*,*) 'lth,lphi= ',lth,lphi
          write(*,*) 'itp=',iitp(q),kind(iitp(q))
          write(titp,"(i0)") iitp(q)
          write(*,*) 'titp=',titp
          result = makedirqq("dat/itp/"//trim(titp))
          result = makedirqq("dat/itp/"//trim(titp)//"/nemato/")
          exit
         endif
         it = it + 1

        enddo
        ppi = acos(-1.d0)
        write(*,*) ' '
        write(*,*) brold(:,71)
        write(*,*) ' '
        write(*,*) bzold(:,0)

        allocate(vecr(-1:lyy))
        allocate(vecth(-1:lth))
        allocate(vecphi(-1:lphi))

        dx = 1.d0/(real(lyy))
        do ix = -1,lyy
           vecr(ix) = real(ix)*dx + dx/2.                                        
        enddo
 
        do ith = -1,lth
           vecth(ith) = real(ith)*2.d0*ppi/real(lth)
        enddo
 
        do iphi = -1,lphi
           vecphi(iphi) = real(iphi)*2.d0*ppi*r0/real(lphi)
        enddo
        write(*,*) 'vechphi',vecphi
 
        allocate(matr(-1:lyy,-1:lth,-1:lphi))
        allocate(matth(-1:lyy,-1:lth,-1:lphi))
        allocate(matphi(-1:lyy,-1:lth,-1:lphi))
        allocate(matb1(-1:lyy,-1:lth,-1:lphi))
        allocate(matb2(-1:lyy,-1:lth,-1:lphi))
        allocate(matb3(-1:lyy,-1:lth,-1:lphi))
 
!        write(cmd,'(a,i0.2,a)') 'mkdir -p dat/itp/', itp,'/nemato'
!        write(*,*) 'CMD',cmd
!        call system(cmd)
        itp = int(iitp(q))
        if (shsh .eq. .false.) then
         tfile = 'dat/itp/'//trim(aint2char(itp))//'/nemato/poinc'
     .        //trim(aint2char(lth))
     .        //'x'//trim(aint2char(lphi))//'.bin'
        else
         tfile = 'dat/itp/'//trim(aint2char(itp))//'/nemato/poinc'
     .         //trim(aint2char(lth))
     .         //'x'//trim(aint2char(lphi))//'_SH_benchmark.bin'
        endif

        open(unit=uft,file=trim(tfile)
     .       ,status='unknown',form='unformatted',action='write')


! Beginning of parallel region
       write(*,*) 'Beginning parallel region'	
!       thr_tot=omp_get_num_threads()
!!!$OMP PARALLEL 
!!!$OMP.   COPYIN(lth,lphi,lyy)
	ieli = 1.d0 * inteli
        write(*,*) 'ieli= ',ieli
!#! parameters for mode_dump
!#! mode not to be dumped
        nmode = 1.d0 * inteli
!#! dump factor
        if (mode_dump) then
          fact = 1.d0/1.d0
        endif

!#! dominant mode dumping
        if (dominant_dump) then
          dom_fact = 1.00
        endif

        write(*,*) 'limit',lphi,lth,lyy
!!!$OMP DO PRIVATE(iphi,ith,ix) SCHEDULE(DYNAMIC) 
        do iphi = -1, lphi
         iphimod = mod(iphi,50) 
         if (iphimod .eq. 0) then 
          write(*,*) 'iphi=',iphi
         endif
!         call flush(6)
         
         aphi = vecphi(iphi) / r0
!         thr_num=omp_get_thread_num()

         do ith = -1,lth
          ath = vecth(ith)

          do ix = -1,lyy
 
           matr(ix,ith,iphi) = 0.0
           matth(ix,ith,iphi) = 0.0
           matphi(ix,ith,iphi) = 0.0
           matb1(ix,ith,iphi) = 0.0
           matb2(ix,ith,iphi) = 0.0
           matb3(ix,ith,iphi) = 0.0
                                                                        
! antitrasformata di fourier:
            do im=0,my
	      m = mm(im) 
              do n = nanf(m),nanf(m)+nzcon(m)-1

               jmn = janz3(m,n,nz,nanf)

!
!!                if (n .ne. 0) cycle
!!                if (n .gt. 0) cycle
!!                 if (real(n) .ne. ieli*m) cycle
!!	      if (n.eq.0) write(*,*) 'n',n                      
!!              if (m .eq. 0) then
!!                 if (n .lt. -23) cycle
!!              elseif (m .eq. 1) then
!!                 if ((n .gt. -1) .or. (n .lt. -20)) cycle
!!!	elseif (m .eq.  2) then
!!!                 if ((n .gt. -13) .or. (n .lt. -19)) cycle
!!!	elseif (m .eq.  3) then
!!!                 if ((n .gt. -21) .or. (n .lt. -27)) cycle
!!!	elseif (m .eq.  4) then
!!!                 if ((n .gt. -29) .or. (n .lt. -35)) cycle
!!              endif
!!                 if (m .eq. 1) then
!!                  if (n .le. -6 .and. n .ge. -12) then
!!                   cycle
!!                  endif
!!                 endif
!!                  if (n .eq. m*ieli) cycle
!	             if (iphi .eq. 0 .and. ith .eq. 0 .and. 
!     .                ix .eq. 0) then
!                      write(*,'(a,2i3,f14.7)') 'm,n',m,n
!                     endif
!
!!#! SINGLE HELICITY
               if (shsh .eq. .true.) then

                if ( n .ne. ieli*m ) then
                 cycle ! prende la stesse elicità di psi
                endif
                if (ix .eq. 0 .and. ith .eq. 0 .and. iphi .eq. 0) then
                 write(*,'(a,4i5)') 'SH: m,n,janz ',ieli,m,n,jmn
                endif

               else

!!               if (m .gt. 1) cycle
!!               if (m .ne. 0) cycle
                 if (ix .eq. 0 .and. ith .eq. 0 .and. iphi .eq. 0) then
                  write(*,'(a,3i5)') 'ALL MODES: m,n,janz ',m,n,jmn
                 endif

               endif
!                nmcheck = nmcheck + 1                                   
!
!
!!#! phases selections (normally zero)
!               if (random_phases) then
!                 if (m .eq. 1 ) then
!!                   if (n.eq.-10 .or. n.eq. -9 .or. n .eq. -8 
!!     .              .or. n .eq. -7 .or. n .eq. -5) then
!                    if (n .eq. -6) then
!                     phase = 2.
!!                    if (n .lt. -3 .and. n .gt. -10) then
!!                     call random_number(rnd)
!!                     phase = rnd*2.*3.141592653589d0
!	             if (iphi .eq. 0 .and. ith .eq. 0 .and. 
!     .                ix .eq. 0) then
!                      write(*,'(a,2i3,f14.7)') 'phase',m,n,phase
!                     endif
!
!                    endif
!                  else 
!                   phase = 0.d0 
!                  endif
!                else
!                 phase = 0.d0
!                endif
!
!               if (ix .eq. 0 .and. ith .eq. 0 .and. iphi .eq. 0) then
!!                write(*,*) 'phase ',phase
!               endif
!
             
               if (m .eq. 0 .and. n .eq. 0) foo = 1.
               if (m .ne. 0 .or.  n .ne. 0) foo = 0.
!                         
!#! Marco, leggo i campi magnetici direttamente dai file output di
!SpeCyl, quindi è necessario aggiungere il fattore 2
                fmns = - 2. * sin (m*ath + n* aphi)
                fmnc =   2. * cos (m*ath + n* aphi)
!
!!#! secondary modes dumping
!              if (mode_dump) then
!               if (m .eq. 0) then
!                if (n .ne. 0) then
!                 abrix = fact*abrix
!!                if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!!                 write(*,'(a,i3,3f9.4)') 'hheck',m,real(n),
!!     .                                -nmode,fact
!!                endif
!                endif
!               endif
!               if (m .eq. 1) then
!                if (real(n) .ne. -nmode) then
!                 abrix = fact*abrix
!!                if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!!                 write(*,'(a,i3,3f9.4,2f13.8)') 'hheck',m,real(n),
!!     .                                -nmode,fact,abrix,abrix/fact
!!                endif
!                endif
!               endif
!               if (m .eq. 2) then
!                if (real(n) .ne. -2*nmode) then
!                 abrix = fact*abrix
!!                if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!!                 write(*,'(a,i3,3f9.4)') 'hheck',m,real(n),
!!     .                                -nmode,fact
!!                endif
!                endif
!               endif
!               if (m .eq. 3) then
!                if (real(n) .ne. -3*nmode) then
!                 abrix = fact*abrix
!!                if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!!                 write(*,'(a,i3,3f9.4)') 'hheck',m,real(n),
!!     .                                -nmode,fact
!!                endif
!                endif
!               endif
!               if (m .eq. 4) then
!                if (real(n) .ne. -4*nmode) then
!                 abrix = fact*abrix
!!                if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!!                 write(*,'(a,i3,3f9.4)') 'hheck',m,real(n),
!!     .                                -nmode,fact
!!                endif
!                endif
!               endif
!              endif
!              
!              if (dominant_dump) then
!               if (m .ne. 0) then
!                if (n .eq. ieli) then
!                 abrix = abrix/dom_fact
!                 if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!                 write(*,'(a,i3,f9.4,2f13.7)') 'dominant_dump',m,
!     .                real(n),abrix,dom_fact
!                  write(*,*) abrix
!                 endif
!                endif
!               endif
!              endif
!
                 abrix = (brold(ix,jmn) + brold(ix+1,jmn) )/2.d0
                 if (ix.eq.-1) then
                  abrix = -(brold(0,jmn) + brold(1,jmn))/2.d0
                 endif
                 if (ix.eq.lx) then 
                  abrix = -(brold(lx-1,jmn) + brold(lx,jmn))/2.d0
                 endif
!             if (n .eq. ieli .and. ith .eq. 0 .and. iphi .eq. 0) then
!               write(*,*) ix,abrix
!              endif
 
                 matb1(ix,ith,iphi)  = matb1(ix,ith,iphi)
     .                + foo * abrix
!     .                + fmnc * real(abrix)
     .                + fmns * abrix
!
!!c..... campo poloidale
!
!
!!#! dominant mode and harmonics dumping
!
!!#! secondary modes dumping
!              if (mode_dump) then
!               if (m .eq. 0) then
!                if (n .ne. 0) then
!                 abt = fact*abt
!                endif
!               endif
!               if (m .eq. 1) then
!                if (real(n) .ne. -nmode) then
!                 abt = fact*abt
!                endif
!               endif
!               if (m .eq. 2) then
!                if (real(n) .ne. -2*nmode) then
!                 abt = fact*abt
!                endif
!               endif
!               if (m .eq. 3) then
!                if (real(n) .ne. -3*nmode) then
!                 abt = fact*abt
!                endif
!               endif
!               if (m .eq. 4) then
!                if (real(n) .ne. -4*nmode) then
!                 abt = fact*abt
!                endif
!               endif
!              endif
!
!              if (dominant_dump) then
!               if (m .ne. 0) then
!                if (n .eq. m*ieli) then
!                 abt =  abt/dom_fact
!                 if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!                  write(*,'(a,i3,f9.4,2f13.8)') 'dominant_dump',m,
!     .                 real(n),abt,dom_fact
!                 endif
!                endif
!               endif
!              endif
!
!!#! Marco, special for borgogno
!!               if (m .eq. 1 .and. real(n) .eq. -10.) then
!!                 abt = fact*abt
!!               endif
!
                 matb2(ix,ith,iphi)  = matb2(ix,ith,iphi)
     .                + foo * btold(ix,jmn) 
     .                + fmnc * btold(ix,jmn)
!     .                + fmns * aimag(btold(ix,jmn))
!                 matb2(ix,ith,iphi)  = matb2(ix,ith,iphi) +
!     .                  foo * abt + fmnc * abt + fonc * abt
!
!!c.....campo toroidale
!
!                 abz = real(bz(ix,jmnm))
!!		 if (m .ge. 1 .and. n .eq. -9*m) then
!!		  abz = -abz
!!		 endif
!
!!#! dominant mode and harmonics dumping
!
!!#! secondary modes dumping
!              if (mode_dump) then
!               if (m .eq. 0) then
!                if (n .ne. 0) then
!                 abz = fact*abz
!                endif
!               endif
!               if (m .eq. 1) then
!                if (real(n) .ne. -nmode) then
!                 abz = fact*abz
!                endif
!               endif
!               if (m .eq. 2) then
!                if (real(n) .ne. -2*nmode) then
!                 abz = fact*abz
!                endif
!               endif
!               if (m .eq. 3) then
!                if (real(n) .ne. -3*nmode) then
!                 abz = fact*abz
!                endif
!               endif
!               if (m .eq. 4) then
!                if (real(n) .ne. -4*nmode) then
!                 abz = fact*abz
!                endif
!               endif
!              endif
!
!              if (dominant_dump) then
!               if (m .ne. 0) then
!               if (n .eq. m*ieli) then
!                abz = abz/dom_fact
!                if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!                write(*,'(a,i3,f9.4,2f13.8)') 'dominant_dump',m,real(n),
!     .                               abz,dom_fact
!                endif
!               endif
!               endif
!              endif
!
!!#! Marco, special for borgogno
!!               if (m .eq. 1 .and. real(n) .eq. -10.) then
!!                 abz = fact*abz
!!               endif
!
!
!             if (n .eq. ieli .and. ix .eq. 0 .and. iphi .eq. 0) then
!              write(*,*) fmnc,fmns,foo
!             endif
                 matb3(ix,ith,iphi)  = matb3(ix,ith,iphi)
     .                + foo * bzold(ix,jmn)
     .                + fmnc * bzold(ix,jmn)
!     .                + fmns * aimag(bzold(ix,jmn))
!              matb3(ix,ith,iphi)  = matb3(ix,ith,iphi) +
!     .                  foo * abz + fmnc * abz + fonc * abz
!
!             if (n .eq. ieli .and. ith .eq. 0 .and. iphi .eq. 0) then
!               write(*,*) ix,abrix,btold(ix,jmn)
!              endif
           enddo
           enddo
          enddo
          enddo
         enddo
!!!$OMP END DO
!!$OMP END PARALLEL
      
        ly2 = int(lyy)
        lth2 = int(lth)
        lphi2 = int(lphi)
        bcond = int(bcond)
        call flush(6)

	write (*,*) 'lxx,lth,lphi,bcond,tor'
	write (uft) ly2,lth2,lphi2,bcond,tor
	write (*,*) 'vecr'
	write (uft) vecr
	write (*,*) 'vecth'
	write (uft) vecth
	write (*,*) 'vecphi'
	write (uft) vecphi
	write (*,*) 'matb1'
	write (uft) matb1
	write (*,*) 'matb2'
	write (uft) matb2
	write (*,*) 'matb3'
	write (uft) matb3


	close(uft)


	deallocate(vecr,vecth,vecphi)
	deallocate(matr,matth,matphi)
	deallocate(matb1,matb2,matb3)
       endif !#! end of the IF (iitp(q>1) > 0)
       enddo !#! end of the cycles on iitp
	
        contains

c     int2char
c     ################################################################
      function aint2char(n) result (chr)

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

      end function aint2char
          function janz3(am,an,anz,ananf)
            
            integer(wp), intent(in) :: am,an
            integer(wp), intent(in) :: anz(:),ananf(:)
            integer(wp) :: aj,mi
!            integer(wp), intent(out) :: janz2

              aj = 0
              if (am .eq. 0 .and. an .eq. 0) then 
               janz3 = 0
              else
               if (am .eq. 0) then 
                aj = anz(1) + an + 1
               else
                do mi=0,am-1 
!              write(*,*) 'inside function janz2',am,an,aj,anz(mi+1)
                 aj = aj + anz(mi+1)
!                 write(*,*) 'ins',anz(mi+1)
                enddo
                aj =  aj + an - ananf(mi+1) + 1
               endif
              endif

              janz3 = aj
          end

        end

