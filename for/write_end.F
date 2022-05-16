
        program write_end
!#! Marco, 11 ottobre 2017
!#! Tentativo di velocizzare la lettura dei file di output di SpeCyl
!leggendoli in Fortran invece che in IDL....
	

	character gran*20,t1*30,t2*20,t3*20,t4*20
	character titp*10,tlth*10,tlphi*10

	integer(8) :: ly,izz
	integer(8) :: max_time,time
        integer :: nfiles(1:9),nfilef(1:9)
        integer :: stats0(13),stats1(13)
	integer :: ios
        real(8) :: r0, zeit
        complex(8), allocatable, dimension(:,:) :: br,bt,bz,vr,vt,vz
        real(8), allocatable, dimension(:) :: t
        integer :: uft=41
        integer :: uft1=42
        integer :: uend=57
	logical :: tor=.false.
        logical :: fasi_sganciate,irrot
        logical :: dominant_dump=.false.
        logical :: random_phases=.false.
        integer :: inlun, nmodes, inteli, d3
        integer :: modenm(0:1)
        character(300) :: path,inputfile,tfile
        character(30) :: ditp, heli, dnth, dnphi
        character(len=555) :: spc2file,nettings_file
        character(300) :: prefix,string,unit_all,unit_end
        namelist /nettings/ nfiles,nfilef
        namelist /spc2/ lx,my,iz,i2d,n2d,m2d
     $   ,fasi_sganciate,irrot

         ninput = iargc() 
         if ( ninput .eq. 1 ) then
            call getarg(1, string); read(string,*) nsec
            write(*,*) 'section number : #',trim(string)
         else
            write(*,*) "usage is: write_end.xxx  NSEC"
         end if
!#! Marco, load default values of the spc2 namelist

!#! preliminary I-O operations
        call getcwd(prefix)
        unit_all = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_all.dat")
        unit_end = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_end.dat")


        spc2file = trim(trim(prefix)//"/for/"//
     .      trim(string)//"/spc2.in")
        open(unit=uft, file=trim(spc2file)
     .     ,status='unknown', iostat=ios)
        if (ios /= 0) then
         write(*,*) 'could not open namelist spc2.in'
         stop
        endif
        read(uft,spc2,iostat=ios)
        if (ios /= 0) then
         write(*,*) 'could not read namelist spc2.in',ios
        else
         close(uft)
        endif

        nettings_file = trim(trim(prefix)//"/for/"//
     .      "/nettings.in")
        open(unit=uft, file=trim(nettings_file)
     .     ,status='unknown', iostat=ios)
        if (ios /= 0) then
         write(*,*) 'could not open namelist nettings.in'
         stop
        endif
        read(uft,nettings,iostat=ios)
        if (ios /= 0) then
         write(*,*) 'could not read namelist nettings.in',ios
        else
         close(uft)
        endif
!#!
        max_time = nfiles(nsec) 
        write(*,*) 'max_time',max_time
!#! default values
        d3 = 1
        r0 = 4.d0
        lx = 100
        my = 4
        iz = 135
        i2d = 0

        write(*,*) 'string',unit_all
        open(unit=uft1, file=trim(unit_all)
     .     ,status='unknown', form='unformatted',iostat=ios)
         if (ios /= 0) then
          write(*,*) 'could not open the inputfile',ios
         endif
         read(uft1,iostat=ios) zeit,ly,izz
         write(*,*) 'ly=',ly, ' izz=',izz
        close(uft1)

!#! allocate arrays
        allocate(br(0:ly,0:izz))
        allocate(vr(0:ly,0:izz))
        allocate(bt(-1:ly,0:izz))
        allocate(vt(-1:ly,0:izz))
        allocate(bz(-1:ly,0:izz))
        allocate(vz(-1:ly,0:izz))
        allocate(t(0:max_time))
        write(*,*) 'allocated'
        
        open(unit=uft, file=trim(unit_all)
     .     ,status='unknown', form='unformatted',iostat=ios)
        open(unit=uend, file=trim(unit_end)
     .     ,status='unknown', form='unformatted',iostat=ios
     .     ,action='write')

!#! ridefinisco max_time per catturare specifici istanti
        do time=0,max_time
         read(uft,iostat=ios) zeit,ly,izz
         if (ios .lt. 0) then
          write(*,*) 'end of file reached',time
          exit
         endif
         t(time) = zeit

!          write(*,*) 'ios',time,ios
         if (time .eq. max_time) then 
          write(uend) zeit,ly,izz
         endif
         call flush(6)
         
         do j=0, izz
          read(uft,iostat=ios) vt(-1,j), vz(-1,j)
     .            , bt(-1,j), bz(-1,j)
          if (time .eq. max_time) then 
           write(uend)  vt(-1,j), vz(-1,j), bt(-1,j), bz(-1,j)
          endif
          if (ios .ne. 0) then
           write(*,*) 'error reading', ios
          endif
          do ix=0, ly
           read(uft,iostat=ios) 
     .              vr(ix,j), vt(ix,j), vz(ix,j)
     .             ,br(ix,j), bt(ix,j), bz(ix,j) 
           if (time .eq. max_time) then 
            write(uend) 
     .              vr(ix,j), vt(ix,j), vz(ix,j)
     .             ,br(ix,j), bt(ix,j), bz(ix,j) 
           endif
           if (ios .ne. 0) then
            write(*,*) 'error reading', ios
           endif
          enddo

!          if (time .gt. 234) then
!           write(*,*) j,bz(time,0,j)
!          endif
!          if (j .eq. 0) then 
!           write(*,*) bz(time,0,j)
!          endif
         enddo
!#!
         if (mod(time,25) .eq. 0) then 
          write(*,*) 'time=',time,t(time)
         endif

        enddo !#! end of dowhile cycle
        write(*,*) 'finished reading! time= ',time
         
        close(uend)
!
!!#! deallocare
        deallocate(t,br,vr,bt,vt,bz,vz)

        end
