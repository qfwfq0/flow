        program exp_poinc
        
        use modpoinc

        use defs
	
        external getarg ! in order to get command line arguments

	character gran*20,t1*30,t2*20,t3*20,t4*20
	character titp*10,tlth*10,tlphi*10

!	integer, parameter :: lth=32
!	integer, parameter :: lphi=64
        character(len=250) :: shot, time
	character(len=250) :: nrr, nth, nphi
        integer(wp) :: nr0,nth0,nphi0,time0,shot0,solo_7
        integer :: nr2,nth2,nphi2
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
!$OMP THREADPRIVATE(aphi,ath)
        real(wp),save :: abrix,abt,abz,fact,dom_fact
!$OMP THREADPRIVATE(abrix,abt,abz,fact,dom_fact)
        real(wp),save :: nmode!, r0
!$OMP THREADPRIVATE(nmode, r0)
        real(wp), allocatable, dimension(:) :: vecr
        real(wp), allocatable, dimension(:)  :: vecth
        real(wp), allocatable, dimension(:)  :: vecphi
!        complex(8),allocatable,dimension(:,:) :: bbr
!        complex(8),allocatable,dimension(:,:) :: bbth
!        complex(8),allocatable,dimension(:,:) :: bbphi
        real(wp), allocatable,dimension(:) :: vecr2,vecth2,vecphi2
        real(wp),allocatable,dimension(:,:,:) :: matr
        real(wp),allocatable,dimension(:,:,:) :: matth
        real(wp),allocatable,dimension(:,:,:) :: matphi
        real(wp),allocatable,dimension(:,:,:) :: matb1
        real(wp),allocatable,dimension(:,:,:) :: matb2
        real(wp),allocatable,dimension(:,:,:) :: matb3
        real(wp),allocatable,dimension(:,:,:) :: matbx
        real(wp),allocatable,dimension(:,:,:) :: matby
        real(wp),allocatable,dimension(:,:,:) :: matbz
        real(wp),allocatable,dimension(:,:,:) :: jac
        real(wp),allocatable,dimension(:,:,:) :: xxc
        real(wp),allocatable,dimension(:,:,:) :: yyc
        real(wp),allocatable,dimension(:,:,:) :: zzc
        real(wp),allocatable,dimension(:,:,:,:) :: car
        integer :: uftt=411
        integer :: uout = 456
	integer :: bcond(6)
        integer(wp) :: section
!	logical :: tor=.false.
	logical :: tor=.true.
        logical,save :: mode_dump
        logical,save :: dominant_dump=.false.
        logical,save :: random_phases=.false.
!$OMP THREADPRIVATE(mode_dump,dominant_dump,random_phases)
	integer, dimension (9):: nfiles
        integer :: inlun, nmodes, inteli, d3
        integer(wp) :: stat
        integer :: modenm(0:1)
        character(300) :: path,inputfile,tfile
        character(30) :: ditp, heli, dnth, dnphi
        character(len=255) :: nmtfile
        character(len=255) :: unit_name,solo7
        integer(wp) :: dumdum

!        namelist /nmt/ itp,inteli,lth,lphi,mode_dump
!#! Marco, load default values of the spc2 namelist
        itp = 0
        inteli = 0
        lth = 32
        lphi = 64
        mode_dump = .false.

        write(*,*) 'tor=',tor
!        shot = 31893
!        time = 150
!        nrr = 101
!        nth = 32
!        nphi = 64
        call getarg(1, string); read(string,*) shot
        call getarg(2, string); read(string,*) time
        call getarg(3, string); read(string,*) nrr
        call getarg(4, string); read(string,*) nth
        call getarg(5, string); read(string,*) nphi
        call getarg(6, string); read(string,*) unit_name
        call getarg(7, string); read(string,*) tfile
        call str2int(time,time0,stat)
!      subroutine str2int(str,int,stat)
!        write(*,*) 'bb',shot,time,nrr,nth,nphi
        call str2int(nrr,nr0,stat)
        call str2int(nth,nth0,stat)
        call str2int(nphi,nphi0,stat)
!        write(*,*) 'bb',shot,time,nrr,nth,nphi,nr0,nth0,nphi0
    
        if ( (trim(shot) .eq. "") .or.
     .     (trim(time) .eq. "")) then
         write(*,*) 'usage: ./pppoincare shot time(in ms) nrr nth nphi
     .inputfilename outfilename'
         stop
        endif  

        bcond(1) = 7
        bcond(2) = 0
        bcond(3) = 1
        bcond(4) = 1
        bcond(5) = 1
        bcond(6) = 1

        d3 = 1
        r0 = 1.995d0
        
        prefix="/0/rfx/ricercatori/ft/specyl/archive/veranda/exp_data/"
        write(*,*) 'lth,lphi= ',trim(prefix)
        write(*,*) 'ccc',trim(shot)
        write(*,*) 'unit=',trim(trim(prefix)//"shot."
     .//trim(shot)//"/"//unit_name)
        open(unit=uftt,file=trim(unit_name),action='read'
     .        ,form='unformatted')

         read(uftt) nr0
         read(uftt) nth0
         read(uftt) nphi0
         allocate(vecr(-1:nr0-2))
         allocate(vecr2(0:nr0-1))
         allocate(vecth(-1:nth0-2))
         allocate(vecth2(0:nth0-1))
         allocate(vecphi(-1:nphi0-2))
         allocate(vecphi2(0:nphi0-1))
         write(*,*) 'nrr2',nr0,nth0,nphi0,shape(vecr)
         read(uftt), vecr
         write(*,*) 'ciao1'
         read(uftt), vecth
         write(*,*) 'ciao2',vecth(-1),vecth(nth0-2)
         write(*,*) 'ciao2',vecth
         read(uftt), vecphi
         write(*,*) 'ciao3',vecphi(-1),vecphi(nphi0-2)
          
!#!     allocate arrays
         allocate(matb1(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         allocate(matb2(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         allocate(matb3(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         read(uftt) matb1
         read(uftt) matb2
         read(uftt) matb3
         allocate(matbx(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         allocate(matby(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         allocate(matbz(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         read(uftt) matbx
         read(uftt) matby
         read(uftt) matbz
         allocate(jac(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         allocate(xxc(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         allocate(yyc(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         allocate(zzc(-1:nr0-2,-1:nth0-2,-1:nphi0-2))
         read(uftt) jac
         read(uftt) xxc
         read(uftt) yyc
         read(uftt) zzc


        close(uftt)

        allocate(car(-1:nr0-2,-1:nth0-2,-1:nphi0-2,3))
        car(:,:,:,1) = xxc
        car(:,:,:,2) = yyc
        car(:,:,:,3) = zzc

!#! open uout
        nr2 = nr0 -2
        nth2 = nth0 -2
        nphi2 = nphi0 -2
!         tfile="/ricercatori/ft/specyl/archive/veranda/
!     .exp_data/shot.31893/poincprova.bin"
!      if (solo_7 .eq. 1) then 
!        tfile=trim(trim(prefix)//"shot."
!     .//trim(shot)//"/"//"poinc_time"//trim(time)
!     .//"_"//trim(nrr)//"x"//trim(nth)
!     .//"x"//trim(nphi)//"_SH_EQ.bin")
!      else
!        tfile=trim(trim(prefix)//"shot."
!     .//trim(shot)//"/"//"poinc_time"//trim(time)
!     .//"_"//trim(nrr)//"x"//trim(nth)
!     .//"x"//trim(nphi)//".bin")
!      endif
        open(unit=uout,file=trim(tfile),status='unknown',
     .  form='unformatted')
	write (*,*) 'lxx,lth,lphi,bcond,tor'
	write (uout) nr2,nth2,nphi2,bcond,tor
	write (*,*) 'vecr'
	write (uout) vecr
	write (*,*) 'vecth'
	write (uout) vecth
	write (*,*) 'vecphi'
	write (uout) vecphi
	write (*,*) 'matb1'
	write (uout) matb1
	write (*,*) 'matb2'
	write (uout) matb2
	write (*,*) 'matb3'
	write (uout) matb3
	write (*,*) 'matbx'
	write (uout) matbx
	write (*,*) 'matby'
	write (uout) matby
	write (*,*) 'matbz'
	write (uout) matbz
	write (*,*) 'jac'
	write (uout) jac
	write (*,*) 'car'
	write (uout) car

	close(uout)


	deallocate(vecr,vecth,vecphi)
	deallocate(matb1,matb2,matb3)
	deallocate(matbx,matby,matbz)
	deallocate(xxc,yyc,zzc,jac)
	


        end
