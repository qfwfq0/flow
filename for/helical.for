        program helical
        
        include 'cyl1.inc.for'

!      inizio ppcon1.for

        character gran*20,t1*30,t2*20,t3*20,t4*20
        character titp*10,tlth*10,tlphi*10,tfile*60

        real(8) :: fmns,fmnc,fons,fonc
        real(8) :: rnd2,eli,pi
        real(8) :: fa,fb,fc
        parameter (lth=8)
        parameter (lphi=8)
!#! attention: n_s_p=a^2
        parameter (n_s_p=400) 
        real(8) :: vecr(-1:lx)
        real(8) :: vecth(-1:lth)
        real(8) :: vecphi(-1:lphi)
        real(8) :: start_point(0:1,0:n_s_p)
        real(8),allocatable,dimension(:) :: helflux0
        real(8),allocatable,dimension(:) :: helflux1
        integer :: uft=41
        integer :: bcond(6)
        integer :: idx0
        logical :: tor=.false.
        	
        integer, dimension (9):: nfiles
        COMMON /COMSET/ nfilei,nfilef,nfiles,nmedi,nmedf,ipag,ilxca

        parameter (nfield=50)
        integer, dimension (100):: nc

        COMMON /COMPOINCARE/ nc

        bcond(1) = 7
        bcond(2) = 0
        bcond(3) = 1
        bcond(4) = 1
        bcond(5) = 1
        bcond(6) = 1
        
        include 'cyl2.inc.for'
        pi = acos(-1d0)
        eli = -7.d0
        dx = 1./lx                                                       
        write(*,'(f17.12)') dx
        DO ix = -1,lx
           vecr(ix) = ix*dx + dx/2.                                        
        ENDDO

        DO ith = -1,lth
           vecth(ith) = ith*2.*PI/lth
        ENDDO

        DO iphi = -1,lphi
           vecphi(iphi) = iphi*2.*PI*RR/lphi
        ENDDO

        allocate(helflux0(-1:lx))
        allocate(helflux1(-1:lx))

! inizia ppcon2.for
        
!       do 500 itp = 1 , maxval(nc)
       do 500 itp = nfilei , maxval(nc)

       t(itp) = zeit

!       write(6,*)'         zeit = ',zeit

! finisce ppcon2.for

        INCLUDE 'cyl3.inc.for'

!  inizia ppcon3.for

! se è un istante prescelto lo processa, altrimenti salta
        iprocess = 0
        do inc=2,size(nc)
           if (nc(inc).ne.0 .and. itp .eq. nc(inc)) iprocess = 1
        enddo
! se è il primo istante va bene anche zero
        if (itp .eq. nc(1)) iprocess = 1
        if (iprocess.eq.0) go to 500

        write(titp,"(I0)") itp
        write(tlth,"(I0)") lth
        write(tlphi,"(I0)") lphi
!	write(*,*) titp

        call system('mkdir -p itp/'//trim(titp))
        call system('mkdir -p itp/' //trim(titp)//'/nemato')


        tfile = 'itp/'//trim(titp)//'/nemato/starting_point.bin'


        open(unit=uft,file=tfile,status='unknown',form='unformatted')

!        DO 501 iphi = -1, lphi
!
!        aphi = vecphi(iphi)/RR
!
!        write(6,*)' iphi = ',iphi,'  phi(iphi) = ',aphi
!
!        ZMAX = -1.E38                                                   
!        ZMIN =  1.E38                                                   
!
!        DO 10 ith = -1,lth
!
!           ath = vecth(ith)
!           athdeg = ath*360./(2.*PI)
!
!!	   write(6,*)' ith = ',ith,'  ath = ',ath,'  ath deg = ',athdeg
!!	   write(6,*)' ith = ',ith
!
        DO 11 ix = -1,lx

           helflux0(ix) = 0.0
           helflux1(ix) = 0.0
                                                                        
        NMCHECK = 0
                                                                        
            DO 30 IM=0,MY
                                                                        
               M = MM(IM) 
              if (m .gt. 1) cycle
!	      if (m .ne. 1 .and. m .ne. 0) cycle ! prende solo i modi con m scelto (solo m minore di 1)

              DO 33 N = NANF(M),NANF(M)+NZCON(M)-1
!	      if (n.eq.0) write(*,*) 'N',N                      
ccc                                                                       
                NMCHECK = NMCHECK + 1                                   
                if (n .ne. eli*m) cycle

!  finisce ppcon3.for

        INCLUDE 'cyl4.inc.for'

!  inizia ppcon4.for
                                                                        
               jmn = janz(m,n)                                          

               if(m .eq. 0 .and. n .eq. 0) foo = 1.
               if(m .ne. 0 .or.  n .ne. 0) foo  =  0.                   
                                                                        
               if(m.eq.0 .and. n.ne.0) then                             
                fons = - 2. * sin (n* aphi)
                fonc =   2. * cos (n* aphi)
               else                                                     
                fons = 0.                                               
                fonc = 0.                                               
               endif                                                    
                         
                if( m .ne. 0 ) then
                 fmns = - 2. * sin (m*ath + n* aphi)
                 fmnc =   2. * cos (m*ath + n* aphi)
                else                                                    
                 fmns =  0.                                          
                 fmnc =  0.                                          
                endif                                                   

!..... helical flux function
                                                        
                ag =  (aaz1+aaz2)/2.
     .		- eli*x2(ix)/RR*(aat1+aat2)/2.

        	 if (ix.eq.-1) ag = 0.d0
!        	 if (ix.eq.lx) ag = -(aimag(br(lx-1,j))+aimag(br(lx,j)))/2.
                 
!                 helflux(ix,ith,iphi) = ag * foo + ag * fonc + ag * fmnc
                 if (m .eq. 0) then
                  helflux0(ix) = (aaz1+aaz2)/2.
     .		- eli*x2(ix)/RR*(aat1+aat2)/2.
                 else if (m .eq. 1) then
                  helflux1(ix) = (aaz1+aaz2)/2.
     .		- eli*x2(ix)/RR*(aat1+aat2)/2.
                 endif
               if (iphi .eq. 0 .and. ith .eq. 0 .and. ix .eq. 1) then
                  write(*,'(a,2i3)') 'before1',m,n
               endif
               if (m .eq. 0 .and. iphi .eq. 0 .and. ith .eq. 0) then
!                  write(*,'(a,i3,2f15.9)') 'before2',ix,x2(ix),
!     .                    helflux0(ix)
               endif
               if (m .eq. 1 .and. iphi .eq. 0 .and. ith .eq. 0) then
!                  write(*,'(a,i3,2f15.9)') 'before3',ix,x2(ix),
!     .                    helflux1(ix,ith,iphi)
               endif
                                                                        
 33	  continue                                                     
 30	continue                                                    

11      continue                                                       
         helflux0(-1) = -helflux0(0)
         helflux1(-1) = - helflux1(0)
         helflux0(lx) = helflux0(lx-1)
         helflux1(lx) = helflux1(lx-1)
!10       continue                                                       
!
!!  ....................................   parametri per la grafica
!                                                                        
!        write(t1,4005) zeit
!4005    format('T = ',f9.1,' $')
!                                                                        
!!        write(6 ,4005) zeit
!
!501    continue
         
         idx0 = 52
         
         fa = helflux0(idx0)
         fb = helflux1(idx0)
         fc = -0.02d0
         write(*,*) fa,fb,fc
         write(*,*) (fc-fa)/fb

         rnd2 = 4.1d0
         do j=0,n_s_p
          call random_number(rnd2)
          start_point(0,j) = rnd2 * 2.d0 * pi
          start_point(1,j) = ( - start_point(0,j) + acos((fc-fa)/fb))
     .                       / real(eli) * RR 
          if (start_point(1,j) .lt. 0.d0) then
           start_point(1,j) = start_point(1,j) + 2.d0*pi
          endif
          write(*,*) j,start_point(0,j),start_point(1,j)
         enddo
         write(uft) n_s_p
         write(uft) x2(idx0)
         write(uft) start_point
!        write (*,*) 'lx,lth,lphi,bcond,tor'
!        write (uft) lx,lth,lphi,bcond,tor
!        write (*,*) 'vecr'
!        write (uft) vecr
!        write (*,*) 'vecth'
!        write (uft) vecth
!        write (*,*) 'vecphi'
!        write (uft) vecphi
!        write (*,*) 'matb1'
!        write (uft) matb1
!        write (*,*) 'matb2'
!        write (uft) matb2
!        write (*,*) 'matb3'
!        write (uft) matb3
!	write (*,*) 'matbx'
!	write (uft) matbx
!	write (*,*) 'matby'
!	write (uft) matby
!	write (*,*) 'matbz'
!	write (uft) matbz
!	write (*,*) 'matjac'
!	write (uft) matjac
!	write (*,*) 'matcar'
!	write (uft) matcar

        close(uft)

500    continue

        deallocate(helflux0,helflux1)
        
       stop
      
        end

      include 'eq_csha.blc.for'
      include 'settings.blc.for'
      include 'pppoincare.blc.for'
