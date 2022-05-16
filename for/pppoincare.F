
        program pppoincare_par
        
        use modpoinc

        use defs
	

	character gran*20,t1*30,t2*20,t3*20,t4*20
	character titp*20,tlth*10,tlphi*10,cmd*100
        character(35):: str

!	integer, parameter :: lth=32
!	integer, parameter :: lphi=64
	integer(4) :: ly2,lth2,lphi2
!;	integer(wp) :: iphimod
	integer(wp),save :: lth,lphi
!$OMP THREADPRIVATE(lth,lphi)
!      OMP variables
	integer,save :: thr_tot,thr_num
!$OMP THREADPRIVATE(thr_tot,thr_num)
	integer :: omp_get_num_threads,omp_get_thread_num
	integer(wp),save :: pm,pjmn
!$OMP THREADPRIVATE(pm,pjmn)
	integer(wp),save :: pim,pn
	integer,save :: iphi,ith,iix
!#! Marco, 20 maggio 2019. Con ifort non bisogna dichiarare l'indice del
!ciclo come threadprivate qui ma farlo al momento della dichiarazione
!del OMP DO. Non ho ancora capito chiaramente perché, ma tant'è.
!!$OMP THREADPRIVATE(iphi,ith,iix)
        integer(8) :: infact,ifact
        real(wp),save :: fmns,fmnc,fons,fonc,foo,ieli,rnd
!$OMP THREADPRIVATE(fmns,fmnc,fons,fonc,foo,ieli,rnd)
        real(wp),save :: aphi,ath
!$OMP THREADPRIVATE(aphi,ath)
        real(wp) :: ppi
        real(wp),save :: abt,abz,fact,dom_fact
!$OMP THREADPRIVATE(abt,abz,fact,dom_fact)
        complex(wp), save :: abrix
!$OMP THREADPRIVATE(abrix)
        real(wp),save :: nmode,r00
!$OMP THREADPRIVATE(nmode, r00)
        real(wp), allocatable, dimension(:) :: vecr
        real(wp), allocatable, dimension(:)  :: vecth
        real(wp), allocatable, dimension(:)  :: vecphi
        real(wp),allocatable,dimension(:,:,:) :: matr
        real(wp),allocatable,dimension(:,:,:) :: matth
        real(wp),allocatable,dimension(:,:,:) :: matphi
        real(wp),allocatable,dimension(:,:,:) :: matb1
        real(wp),allocatable,dimension(:,:,:) :: matb2
        real(wp),allocatable,dimension(:,:,:) :: matb3
!        integer :: uft=41
	integer(4) :: bcond(6)
        integer(wp) :: section
        logical :: result
	logical :: tor=.false.
        logical,save :: mode_dump
!$OMP THREADPRIVATE(mode_dump)
        logical,save :: shsh
        logical,save :: dominant_dump=.false.
        logical,save :: random_phases=.false.
!$OMP THREADPRIVATE(dominant_dump,random_phases)
	integer, dimension (9):: nfiles
        integer :: inlun, nmodes,  d3,len_itp
        integer(wp),save :: inteli
!$OMP THREADPRIVATE(inteli)
        integer :: modenm(0:1)
        integer, dimension(0:10) :: iitp
        character(300) :: path,inputfile,tfile
        character(300) :: prefix2,prefix0
        character(30) :: ditp, heli, dnth, dnphi
        character(len=255) :: nmtfile

        namelist /nmt/ iitp,inteli,lth,lphi,shsh,mode_dump,fact
!#! Marco, load default values of the spc2 namelist
        iitp(0:10) = 0
        inteli = 0
        lth = 32
        lphi = 64
        mode_dump = .false.
        shsh = .false.
        fact = 1.d0

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
!        call getcwd(prefix2)
!        pref = index(prefix2,'/ricercatori/')
!        write(*,*) 'AAA',prefix2
!        write(*,*) 'bbb',pref
!        prefix = trim("/0/rfx"//prefix2(pref:len(prefix2)))
        call getcwd(prefix)
        write(*,*) 'prefix= ',trim(prefix)
!        prefix='/0/rfx/mnt/N61a/ricercatori_exp4/ft/specyl/
!     .archive2/kryzhanovskyy/tokamak/plasmoidi/
!     .s1e8p100_alpha0.625_
!     .lx256_nmodes32_dt1d-4_pertvr_fasi_sganciate_nuovo_
!     .quat_1_2_3_parts'
!        unit_all = trim(trim(prefix)//"/dat/"
!     .           //"specyl_"//trim(string)//"_all.dat")
  
  !#! Marco, load default values of the spc2 namelist
!        call default_spc2 
   
        spc2file = trim(trim(prefix)//"/for/1/spc2.in")
        spc3file = trim(trim(prefix)//"/for/1/spc3.in")
  
        write(*,*) 'spc2in',trim(spc2file)
        write(*,*) 'spc3in',trim(spc3file)
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
        r00 = 4.d0
        write(*,*) 'SH=',shsh
        write(*,*) 'inteli=',inteli
        write(*,*) 'times to be computed=',iitp
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

           read(ibein,iostat=ierr) zeit,ly,izz

!           if (abs(mod(zeit,10000.)) .lt. 1.d0) then
!            write(*,'(a,f14.7,a,i4,a,i4)') 
!     .         'zeit=',zeit,' ly=',ly,'  izz=',izz
!           endif
            write(*,'(a,f14.7,a,i4,a,i4)') 
     .         'zeit=',zeit,' it=',it,'  err=',ierr
           if (ierr .lt. 0) then !end of file reached
            close (ibein)
            it = it - 1
            section = section + 1
            write(*,*) 'ho aperto lo specyl_?.dat',section
            unit_all = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(aint2char(section))//"_all.dat")
            open(unit=ibein,file=trim(unit_all),action='read'
     .        ,form='unformatted')
            read(ibein,iostat=ierr) zeit,ly,iz
!            if (mod(zeit,100.) .eq. 0.d0) then
!             write(*,*) 'file successivo, section=',section
!             write(*,'(a,f9.5,a,i4,a,i4)') 
!     .         'zeit=',zeit,' ly=',ly,'  izz=',izz
!            endif
           endif

         call flush(6)
         if( izz .ne. iz ) then
          write(6,*) 'numero diverso di modi elaborati',iz,izz
          stop
         endif
         call flush(6)
       
         do j=0,izz
          read(ibein) vt(-1,j), vz(-1,j), bt(-1,j), bz(-1,j)
           do iix=0,ly                                                   
            read(ibein) vr(iix,j), vt(iix,j), vz(iix,j),
     .                 br(iix,j), bt(iix,j), bz(iix,j) 
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
        write(*,*) ly,lth,lphi
!        write(*,*) br(:,11)
!        write(*,*) ' '

        allocate(vecr(-1:ly))
        allocate(vecth(-1:lth))
        allocate(vecphi(-1:lphi))

        dx = 1.d0/(real(ly))
        do iix = -1,ly
           vecr(iix) = real(iix)*dx + dx/2.                                        
        enddo
 
        do ith = -1,lth
           vecth(ith) = real(ith)*2.d0*ppi/real(lth)
        enddo
 
        do iphi = -1,lphi
           vecphi(iphi) = real(iphi)*2.d0*ppi*r00/real(lphi)
        enddo
 
        allocate(matr(-1:ly,-1:lth,-1:lphi))
        allocate(matth(-1:ly,-1:lth,-1:lphi))
        allocate(matphi(-1:ly,-1:lth,-1:lphi))
        allocate(matb1(-1:ly,-1:lth,-1:lphi))
        allocate(matb2(-1:ly,-1:lth,-1:lphi))
        allocate(matb3(-1:ly,-1:lth,-1:lphi))
 
!#! dump factor
        if (mode_dump) then
          infact = int(1./fact)
          ifact = int(fact)
        endif

        itp = int(iitp(q))
        write(cmd,'(a)') 'mkdir -p dat/itp/'
     .//trim(aint2char(itp))//'/nemato'
        call system(cmd)
        if (shsh .eq. .false.) then
         if (mode_dump) then
          write(str,"(f0.3)") fact
          write(*,*) 'aaa',fact,trim(str)
          tfile = 'dat/itp/'//trim(aint2char(itp))//'/nemato/poinc'
     .        //trim(aint2char(lth))
     .        //'x'//trim(aint2char(lphi))//'_sec'//trim(str)//'.bin'
          write(*,*) 'Sto dumpando i modi di un fattore: ',fact
         else 
          tfile = 'dat/itp/'//trim(aint2char(itp))//'/nemato/poinc'
     .        //trim(aint2char(lth))
     .        //'x'//trim(aint2char(lphi))//'.bin'
         endif
        else
         tfile = 'dat/itp/'//trim(aint2char(itp))//'/nemato/poinc'
     .         //trim(aint2char(lth))
     .         //'x'//trim(aint2char(lphi))//'_SH.bin'
        endif

        write(*,*) 'CMD',tfile
        open(unit=uft,file=trim(tfile)
     .       ,status='unknown',form='unformatted',action='write')

! Beginning of parallel region
       write(*,*) 'Beginning parallel region',ieli
       ieli = 1.d0 * inteli
!#! parameters for mode_dump
!#! mode not to be dumped
       nmode = 1.d0 * inteli

!#! dominant mode dumping
       if (dominant_dump) then
        dom_fact = 1000.00
       endif
!       if (thr_num .eq. 0) then
!        write(*,*) 'limit',lphi,lth,ly,thr_num,thr_tot
!       endif

!$OMP PARALLEL NUM_THREADS(8)
!$OMP. COPYIN(lth,lphi,inteli,ieli,r00,mode_dump)
       thr_tot=omp_get_num_threads()
       thr_num=omp_get_thread_num()
!$OMP DO PRIVATE(iphi,ith,iix,pim,pm,pn,aphi,ath,mode_dump)
!$OMP. SCHEDULE(DYNAMIC)
       do iphi = -1, lphi
!        if (thr_num .eq. 2) then 
!        endif
!        if (mod(iphi,50) .eq. 0 .and. thr_num .eq. 0) then 
        if (mod(iphi,50) .eq. 0) then 
         write(*,*) 'iphi=',iphi,'/',lphi
        endif
        aphi = vecphi(iphi) / r00

        do ith = -1,lth
         ath = vecth(ith)

         do iix = -1,ly

          matr(iix,ith,iphi) = 0.0
          matth(iix,ith,iphi) = 0.0
          matphi(iix,ith,iphi) = 0.0
          matb1(iix,ith,iphi) = 0.0
          matb2(iix,ith,iphi) = 0.0
          matb3(iix,ith,iphi) = 0.0
                                                                        
! antitrasformata di fourier:
          do pim=0,my
           pm = mm(pim) 
           do pn = nanf(pm),nanf(pm)+nzcon(pm)-1
            pjmn = janz3(pm,pn,nz,nanf)

!                if (pm .gt. 2) cycle
!                if (pn .ne. 0) cycle
!!                if (pn .gt. 0) cycle
!!                 if (real(pn) .ne. ieli*pm) cycle
!!	      if (pn.eq.0) write(*,*) 'pn',pn                      
!!              if (pm .eq. 0) then
!!                 if (pn .lt. -23) cycle
!!              elseif (pm .eq. 1) then
!!                 if ((pn .gt. -1) .or. (pn .lt. -20)) cycle
!!!	elseif (pm .eq.  2) then
!!!                 if ((pn .gt. -13) .or. (pn .lt. -19)) cycle
!!!	elseif (pm .eq.  3) then
!!!                 if ((pn .gt. -21) .or. (pn .lt. -27)) cycle
!!!	elseif (pm .eq.  4) then
!!!                 if ((pn .gt. -29) .or. (pn .lt. -35)) cycle
!!              endif
!!                 if (pm .eq. 1) then
!!                  if (pn .le. -6 .and. pn .ge. -12) then
!!                   cycle
!!                  endif
!!                 endif
!!                  if (pn .eq. pm*ieli) cycle
!	             if (iphi .eq. 0 .and. ith .eq. 0 .and. 
!     .                iix .eq. 0) then
!                      write(*,'(a,2i3,f14.7)') 'pm,pn',pm,pn
!                     endif
!
!!#! SINGLE HELICITY
            if (shsh .eq. .true.) then
             if ( pn .ne. ieli*pm ) then
              cycle ! prende la stesse elicità di psi
             endif
             if (iix .eq. 0 .and. ith .eq. 0 .and. iphi .eq. 0) then
!     .           .and. thr_num .eq. 0) then
              write(*,'(a,4i5)') 'SH: pm,pn,janz ',inteli,pm,pn,pjmn
             endif
            else
             if (iix .eq. 0 .and. ith .eq. 0 .and. iphi .eq. 0) then
              write(*,'(a,3i5)') 'ALL MODES: pm,pn,janz ',pm,pn,pjmn
             endif
            endif
!                nmcheck = nmcheck + 1                                   
!
!
!!#! phases selections (normally zero)
!               if (random_phases) then
!                 if (pm .eq. 1 ) then
!!                   if (pn.eq.-10 .or. pn.eq. -9 .or. pn .eq. -8 
!!     .              .or. pn .eq. -7 .or. pn .eq. -5) then
!                    if (pn .eq. -6) then
!                     phase = 2.
!!                    if (pn .lt. -3 .and. pn .gt. -10) then
!!                     call random_number(rnd)
!!                     phase = rnd*2.*3.141592653589d0
!	             if (iphi .eq. 0 .and. ith .eq. 0 .and. 
!     .                iix .eq. 0) then
!                      write(*,'(a,2i3,f14.7)') 'phase',pm,pn,phase
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
!               if (iix .eq. 0 .and. ith .eq. 0 .and. iphi .eq. 0) then
!!                write(*,*) 'phase ',phase
!               endif
!
             
            if (pm .eq. 0 .and. pn .eq. 0) then
             foo = 1.d0
             fmns = 0.d0
             fmnc = 0.d0
            else
             foo = 0.d0
             fmnc =   2. * cos (pm*ath + pn* aphi)
             fmns = - 2. * sin (pm*ath + pn* aphi)
            endif
!            if (pim .eq. 0 .and. pn .eq. 0) then
!             write(*,*) 'aaaa',fmnc,fmns
!            endif
!                         
!#! Marco, leggo i campi magnetici direttamente dai file output di
!SpeCyl, quindi è necessario aggiungere il fattore 2
!
            abrix = (br(iix,pjmn) + br(iix+1,pjmn) )/2.d0
            if (iix.eq.-1) then
             abrix = -(br(0,pjmn) + br(1,pjmn))/2.d0
            endif
            if (iix.eq.lx) then 
             abrix = -(br(lx-1,pjmn) + br(lx,pjmn))/2.d0
            endif
 
            if (mode_dump) then
             if (pn .ne. pm*ieli) then
              abrix = abrix * fact
              if (iix.eq.0 .and. ith.eq.0 .and. iphi.eq.0) then
               write(*,'(a,3i5)') 'DUMPED MODES: pm,pn,janz ',pm,pn,pjmn
              endif
             endif
            endif

            matb1(iix,ith,iphi)  = matb1(iix,ith,iphi)
     .                + foo * real(abrix) 
     .                + fmnc * real(abrix)
     .                + fmns * aimag(abrix)
!
!!c..... campo poloidale
!
!
!!#! dominant mode and harmonics dumping
                 if (mode_dump) then
                  if (pn .ne. pm*ieli) then
                   bt(iix,pjmn) = bt(iix,pjmn) * fact
                  endif
                 endif
!
                 matb2(iix,ith,iphi)  = matb2(iix,ith,iphi)
     .                + foo * real(bt(iix,pjmn)) 
     .                + fmnc * real(bt(iix,pjmn))
     .                + fmns * aimag(bt(iix,pjmn))
!                 matb2(iix,ith,iphi)  = matb2(iix,ith,iphi) +
!     .                  foo * abt + fmnc * abt + fonc * abt
!
!!c.....campo toroidale
!!#! dominant mode and harmonics dumping
!
                 if (mode_dump) then
                  if (pn .ne. pm*ieli) then
                   bz(iix,pjmn) = bz(iix,pjmn) * fact
                  endif
                 endif

                 matb3(iix,ith,iphi)  = matb3(iix,ith,iphi)
     .                + foo * real(bz(iix,pjmn))
     .                + fmnc * real(bz(iix,pjmn))
     .                + fmns * aimag(bz(iix,pjmn))
!              matb3(iix,ith,iphi)  = matb3(iix,ith,iphi) +
!     .                  foo * abz + fmnc * abz + fonc * abz
!
           enddo
           enddo
          enddo
          enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      
        ly2 = int(ly)
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

