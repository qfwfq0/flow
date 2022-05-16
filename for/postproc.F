
        program postproc
!#! Marco, 11 ottobre 2017
!#! Tentativo di velocizzare la lettura dei file di output di SpeCyl
!leggendoli in Fortran invece che in IDL....
	

	character gran*20,t1*30,t2*20,t3*20,t4*20
	character titp*10,tlth*10,tlphi*10

	integer(8) :: ly,izz
	integer(8) :: max_time,time
	integer(8) :: time_dat,time_idl
        integer :: nfiles(1:9),nfilef(1:9)
        integer :: stats0(13),stats1(13)
	integer :: ios
        real(8) :: r0, zeit
        complex(8), allocatable, dimension(:,:,:) :: br,bt,bz,vr,vt,vz
        complex(8), allocatable, dimension(:,:,:) :: br0,bt0,bz0,vr0
        complex(8), allocatable, dimension(:,:,:) :: vt0,vz0
        real(8), allocatable, dimension(:) :: t,t0
        integer :: uft=41
        integer :: uft1=42
        integer :: uidl=57
        integer :: uend=157
	logical :: tor=.false.
        logical :: fasi_sganciate,irrot
        logical :: dominant_dump=.false.
        logical :: random_phases=.false.
        integer :: inlun, nmodes, inteli, d3
        integer :: modenm(0:1)
        character(300) :: path,inputfile,tfile
        character(30) :: ditp, heli, dnth, dnphi
        character(len=555) :: spc2file,nettings_file
        character(300) :: prefix,string,unit_all,unit_idl,unit_end
        namelist /nettings/ nfiles,nfilef
        namelist /spc2/ lx,my,iz,i2d,n2d,m2d
     $   ,fasi_sganciate,irrot

         ninput = iargc() 
         if ( ninput .eq. 1 ) then
            call getarg(1, string); read(string,*) nsec
            write(*,*) 'section number : #',trim(string)
         else
            write(*,*) "usage is: "
         end if
!#! Marco, load default values of the spc2 namelist

!#! preliminary I-O operations
        call getcwd(prefix)
        unit_all = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_all.dat")
        unit_idl = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_idl.dat")
        unit_end = trim(trim(prefix)//"/dat/"
     .           //"specyl_"//trim(string)//"_end.dat")


!        time_dat = call system('stat -c%Y'//unit_all//"'")
!        time_idl = call system('stat -c%Y'//unit_idl//"'")
!        write(*,*) time_dat, time_idl
!        stop

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
        allocate(br(0:max_time,0:ly,0:izz))
        allocate(vr(0:max_time,0:ly,0:izz))
        allocate(bt(0:max_time,-1:ly,0:izz))
        allocate(vt(0:max_time,-1:ly,0:izz))
        allocate(bz(0:max_time,-1:ly,0:izz))
        allocate(vz(0:max_time,-1:ly,0:izz))
        allocate(t(0:max_time))
        write(*,*) 'allocated'
        
        open(unit=uft, file=trim(unit_all)
     .     ,status='unknown', form='unformatted',iostat=ios)

        time = 0
!        do while (1 .eq. 1)
        do while (time .le. max_time)
         read(uft,iostat=ios) zeit,ly,izz
         if (ios .lt. 0) then
          write(*,*) 'end of file reached'
          exit
         endif
         t(time) = zeit

!          write(*,*) 'ios',time,ios
         call flush(6)
         
         do j=0, izz
          read(uft,iostat=ios) vt(time,-1,j), vz(time,-1,j)
     .            , bt(time,-1,j), bz(time,-1,j)
          if (ios .ne. 0) then
           write(*,*) 'error reading', ios
          endif
          do ix=0, ly
           read(uft,iostat=ios) 
     .              vr(time,ix,j), vt(time,ix,j), vz(time,ix,j)
     .             ,br(time,ix,j), bt(time,ix,j), bz(time,ix,j) 
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
         time = time + 1 

        enddo !#! end of dowhile cycle
        write(*,*) 'finished reading! time= ',time
         
!#! resize vectors
!#! allocate arrays
        allocate(br0(0:time,0:ly,0:izz))
        allocate(vr0(0:time,0:ly,0:izz))
        allocate(bt0(0:time,-1:ly,0:izz))
        allocate(vt0(0:time,-1:ly,0:izz))
        allocate(bz0(0:time,-1:ly,0:izz))
        allocate(vz0(0:time,-1:ly,0:izz))
        allocate(t0(0:time))
        write(*,*) 'allocated second time',time

        t0 = t(0:time)
        br0 = br(0:time,0:ly,0:izz)
        bt0 = bt(0:time,-1:ly,0:izz)
        bz0 = bz(0:time,-1:ly,0:izz)
        vr0 = vr(0:time,0:ly,0:izz)
        vt0 = vt(0:time,-1:ly,0:izz)
        vz0 = vz(0:time,-1:ly,0:izz)

        write(*,*) 'resized vectors up to time step=',time
        deallocate(t,br,vr,bt,vt,bz,vz)
!#! salvare!
        open(unit=uidl, file=trim(unit_idl)
     .     ,status='unknown', form='unformatted',iostat=ios
     .     ,action='write')
        write(*,*) 'Opened the file for saving'
         write(uidl) ly
         write(uidl) izz
         write(uidl) max_time
         write(uidl) time
         write(uidl) t0
         write(uidl) br0
         write(uidl) bt0
         write(uidl) bz0
         write(uidl) vr0
         write(uidl) vt0
         write(uidl) vz0

        close(uidl)
        write(*,*) 'Closed the file for saving'

!#! deallocare
        deallocate(t0,br0,vr0,bt0,vt0,bz0,vz0)

        end
