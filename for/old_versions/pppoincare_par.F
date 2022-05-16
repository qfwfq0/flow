
        program pppoincare_par
	

	character gran*20,t1*30,t2*20,t3*20,t4*20
	character titp*10,tlth*10,tlphi*10

!	integer, parameter :: lth=32
!	integer, parameter :: lphi=64
	integer :: lth
	integer :: lphi
!$OMP THREADPRIVATE(lth,lphi)
!      OMP variables
	integer :: thr_tot,thr_num
!$OMP THREADPRIVATE(thr_tot,thr_num)
	integer :: omp_get_num_threads,omp_get_thread_num
	integer :: NMCHECK,im,n,m
!$OMP THREADPRIVATE(NMCHECK,im,n,m)
	integer :: iphi,ith,ix
!$OMP THREADPRIVATE(iphi,ith,ix)
        real(8) fmns,fmnc,fons,fonc,foo,ieli,rnd
!$OMP THREADPRIVATE(fmns,fmnc,fons,fonc,foo,ieli,rnd)
        real(8) :: aphi,ath,phase
!$OMP THREADPRIVATE(aphi,ath,phase)
        real(8) :: ZMAX,ZMIN
!$OMP THREADPRIVATE(ZMAX,ZMIN)
        real(8) :: abrix,abt,abz,fact,dom_fact
!$OMP THREADPRIVATE(abrix,abt,abz,fact,dom_fact)
	integer :: jmnm
!$OMP THREADPRIVATE(jmnm)
        real(8) :: nmode, r0
!$OMP THREADPRIVATE(nmode, r0)
        real(8), allocatable, dimension(:) :: vecr
        real(8), allocatable, dimension(:)  :: vecth
        real(8), allocatable, dimension(:)  :: vecphi
        complex(8),allocatable,dimension(:,:) :: bbr
        complex(8),allocatable,dimension(:,:) :: bbth
        complex(8),allocatable,dimension(:,:) :: bbphi
        real(8),allocatable,dimension(:,:,:) :: matr
        real(8),allocatable,dimension(:,:,:) :: matth
        real(8),allocatable,dimension(:,:,:) :: matphi
        real(8),allocatable,dimension(:,:,:) :: matb1
        real(8),allocatable,dimension(:,:,:) :: matb2
        real(8),allocatable,dimension(:,:,:) :: matb3
        integer :: uft=41
	integer :: bcond(6)
	logical :: tor=.false.
        logical :: mode_dump
        logical :: dominant_dump=.false.
        logical :: random_phases=.false.
!$OMP THREADPRIVATE(mode_dump,dominant_dump,random_phases)
	integer, dimension (9):: nfiles
        integer :: inlun, nmodes, inteli, d3
        integer :: modenm(0:1)
        character(300) :: path,string,inputfile,tfile
        character(30) :: ditp, heli, dnth, dnphi
        character(len=255) :: nmtfile

        namelist /nmt/ itp,inteli,lth,lphi,mode_dump
!#! Marco, load default values of the spc2 namelist
        itp = 0
        inteli = 0
        lth = 32
        lphi = 64
        mode_dump = .false.

!#! preliminary I-O operations
        call getcwd(prefix)
        write(*,*) 'prefix',prefix
        unit_all = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_all.dat")

        nmtfile = trim(trim(prefix)//"/for/"//
     .      trim(string)//"/spc2.in")

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
        read( heli, '(i10)' ) inteli
        read( dnth, '(i10)' ) lth
        read( dnphi, '(i10)' ) lphi
        write(*,*) 'lth,lphi= ',lth,lphi

        exit
        call system('mkdir -p ../dat/itp/'//trim(ditp)//'/nemato')
!#! Marco, opening and reading the file containing the magnetic field
!profiles
        inputfile = '../dat/itp/'//trim(ditp)//'/bcyl.bin'
        inlun = 254
        open(inlun,file=inputfile,status='old'
     .       ,form='unformatted')
         read(inlun) lxx
         read(inlun) nmodes
        
         write(*,*) 'lxx, nmodes', lxx, nmodes
         allocate(vecr(-1:lxx))
         allocate(vecth(-1:lth))
         allocate(vecphi(-1:lphi))
         allocate(bbr(0:lxx,nmodes))
         allocate(bbth(0:lxx,nmodes))
         allocate(bbphi(0:lxx,nmodes))
 
         read(inlun) bbr 
         read(inlun) bbth
         read(inlun) bbphi
        close(inlun)

        dx = 1.d0/(real(lxx)-1.d0)                                                      
        do ix = -1,lxx
           vecr(ix) = real(ix)*dx + dx/2.                                        
        enddo
 
        do ith = -1,lth
           vecth(ith) = real(ith)*2.d0*PI/real(lth)
        enddo
 
        do iphi = -1,lphi
           vecphi(iphi) = real(iphi)*2.d0*PI*r0/real(lphi)
        enddo
 
        allocate(matr(-1:lxx,-1:lth,-1:lphi))
        allocate(matth(-1:lxx,-1:lth,-1:lphi))
        allocate(matphi(-1:lxx,-1:lth,-1:lphi))
        allocate(matb1(-1:lxx,-1:lth,-1:lphi))
        allocate(matb2(-1:lxx,-1:lth,-1:lphi))
        allocate(matb3(-1:lxx,-1:lth,-1:lphi))
 
        tfile = '../dat/itp/'//trim(ditp)//'/nemato/poinc'
     .       //trim(dnth)//'x'//trim(dnphi)//'.bin'
        write(*,*) trim(tfile)

        open(unit=uft,file=trim(tfile)
     .       ,status='unknown',form='unformatted')


! Beginning of parallel region
       write(*,*) 'Beginning parallel region'	
       thr_tot=omp_get_num_threads()
!$OMP PARALLEL NUM_THREADS(8)
	ieli = -1.d0 * inteli
!#! parameters for mode_dump
!#! mode not to be dumped
        nmode = 1.d0 * inteli
!#! dump factor
        if (mode_dump) then
          fact = 1.d0/2.25d0
        endif
!#! Marco, special for borgogno
!        fact = 0.8d0

!#! dominant mode dumping
        if (dominant_dump) then
          dom_fact = 1000.00
        endif

!$OMP DO SCHEDULE(DYNAMIC) 
        do iphi = -1, lphi
         if (mod(iphi,10) .eq. 0) write(*,*) 'iphi=',iphi
         call flush(6)
         aphi = vecphi(iphi)/r0

         thr_num=omp_get_thread_num()
!	 write(*,*) 'iphi,thr_num',iphi,aphi,thr_num

         ZMAX = -1.E38                                                   
         ZMIN =  1.E38                                                   

         do ith = -1,lth
          ath = vecth(ith)

          do ix = -1,lxx
 
           matr(ix,ith,iphi) = 0.0
           matth(ix,ith,iphi) = 0.0
           matphi(ix,ith,iphi) = 0.0
           matb1(ix,ith,iphi) = 0.0
           matb2(ix,ith,iphi) = 0.0
           matb3(ix,ith,iphi) = 0.0
                                                                        
           NMCHECK = 0
                                                                        
! antitrasformata di fourier:
            do im=0,nmodes-1
             modenm = mnum(im,nz,nanf,my,d3)
             m = modenm(0)
             n = modenm(1)

!                                                                        
!!               if (m .gt. 1) cycle
!!               if (m .ne. 0) cycle
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
!!		if ( n .ne. ieli*m ) cycle ! prende la stesse elicità di psi
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
               if(m .eq. 0 .and. n .eq. 0) foo = 1.
               if(m .ne. 0 .or.  n .ne. 0) foo  =  0.                   
!                                                                        
!               if(m.eq.0 .and. n.ne.0) then                             
!                fons = - 2. * sin (n* aphi + phase)
!                fonc =   2. * cos (n* aphi + phase)
!!                fons = -  sin (n* aphi)
!!                fonc =    cos (n* aphi)
!               else                                                     
!                fons = 0.                                               
!                fonc = 0.                                               
!               endif                                                    
!                         
!               if( m .ne. 0 ) then
!#! Marco, i campi magnetici provengono da isavemodeprofiles.pro, in cui
!c'è già il fattore 2. dell'antitrasformata di Fourier.
                fmns = - 1. * sin (m*ath + n* aphi + phase)
                fmnc =   1. * cos (m*ath + n* aphi + phase)
!!                fmns = - sin (m*ath + n* aphi)
!!                fmnc =   cos (m*ath + n* aphi)
!               else                                                    
!                fmns =  0.                                          
!                fmnc =  0.                                          
!               endif                                                   
!
!	       if (iphi .eq. 0 .and. ith .eq. 0 .and. 
!     .             ix .eq. 0) then
!	        write(*,*) 'official m,n,jmnm',m,n,jmnm,thr_num
!	       endif
!
!!..... campo radiale
!      
!
!!               if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!!                  write(*,'(a,2i3,f9.4)') 'before1',thr_num,m,real(n)
!!                  write(*,'(a,i3,2f13.8)') 'before2',jmnm,vecr(ix),abrix
!!!                 write(*,'(i2,a,i3,2f9.4,f13.8)') thr_num,' before',m,real(n),
!!!     .                                fact,abrix
!!	       endif
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
!                if (n .eq. m*ieli) then
!                 abrix = abrix/dom_fact
!                 if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!                 write(*,'(a,i3,f9.4,2f13.7)') 'dominant_dump',m,
!     .                real(n),abrix,dom_fact
!                 endif
!                endif
!               endif
!              endif
!
!!#! Marco, special for borgogno
!!               if (m .eq. 1 .and. real(n) .eq. -10.) then
!!                 write(*,*) 'm1,n10 ',m,n,fact
!!                 abrix = fact*abrix
!!               endif
!
!                if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
!!                 write(*,'(a,i3,f9.4,f13.8)') 'final',m,real(n),
!!     .                                abrix
!!                 write(*,*) " "
!                endif
                 matb1(ix,ith,iphi)  = matb1(ix,ith,iphi)
     .                + foo * real(bbr(ix,im)) 
     .                + fmnc * real(bbr(ix,im))
     .                + fmns * aimag(bbr(ix,im))
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
     .                + foo * real(bbth(ix,im)) 
     .                + fmnc * real(bbth(ix,im))
     .                + fmns * aimag(bbth(ix,im))
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
                 matb3(ix,ith,iphi)  = matb3(ix,ith,iphi)
     .                + foo * real(bbphi(ix,im))
     .                + fmnc * real(bbphi(ix,im))
     .                + fmns * aimag(bbphi(ix,im))
!              matb3(ix,ith,iphi)  = matb3(ix,ith,iphi) +
!     .                  foo * abz + fmnc * abz + fonc * abz
!
           enddo
          enddo
          enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL

                                                                        
        call flush(6)
	write (*,*) 'lxx,lth,lphi,bcond,tor'
	write (uft) lxx,lth,lphi,bcond,tor
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
	
!!       stop
        contains
!
!        subroutine prova_janzm(pm,pn,janzm)
!!#! Routine to obtain the right index of the modes
!
!	!Global variables
!	integer, INTENT(IN) :: pm,pn
!
!	!Local variables
!	integer, dimension(:,:), allocatable :: tmp_janzm
!	integer, INTENT(OUT) :: janzm
!	integer :: minm,minn,maxn
!        janzm = 0
!        if (pm .eq. 0 .and. pn .eq. 0) then
!         janzm = 0
!         return
!        endif
!        do im = 0, pm -1 
!         janzm = janzm + nz(im)
!        enddo
!        janzm = janzm + pn - nanf(pm) + 1  
!        end subroutine prova_janzm
!
!;#; this function converts janz to a diad (m,n)
         function mnum(j,nz,nanf,my,d3)
          
           implicit none
           
           integer, intent(in) :: d3, j,my
           integer(8), intent(in) :: nz(0:my)
           integer, intent(in) :: nanf(0:my)
           integer, intent(out) :: mnum(0:1)
           integer :: sizea, mm,q,nn
           integer :: cnz(0:my)
            
           mnum(0) = 0
           mnum(1) = 1
           if (d3 .eq. 1) then 
!          ;#; compute the cumulated nz array, cnz
              if (j .eq. 0) then 
               mnum(0) = 0 ; mnum(1) = 0
              else
              cnz = nz
              sizea=size(nz)
              do q=0, sizea-1
               if (q .eq. 0) then
                cnz(q) = nz(0) - 1
               else
                if (q .eq. 1) then
                 cnz(q) = cnz(q-1) + nz(q) 
                else
                 cnz(q) = cnz(q-1) + nz(q) 
                endif
               endif
              enddo
              do q=0, sizea-1
               if (j-1 .le. cnz(q)) then
                mnum(0) = q
                exit
               endif
              enddo
              if (mnum(0) .eq. 0) then
               mnum(1) = nanf(mnum(0)) + j - 1
              else
               mnum(1) = nanf(mnum(0)) +
     .          j - 1 - cnz(mnum(0)-1) - 1
              endif
              endif
!              write(*,*) j,mnum
            else
!          ;#; 2D case
              if (j .eq. 0) then
              mnum(0) = 0 ; mnum(1) = 0
              else
               mnum(0) = j ; mnum(1) = nanf(j)
              endif
            endif
          end function 
        end

