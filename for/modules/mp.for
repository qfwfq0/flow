      module mp

      use defs
      use math

!#! Marco, subroutines and functions
       contains
        subroutine apply_irrot_bc(mp_percr,amp_ar,m_mpr,n_mpr
     .                           ,ph_mpr,bbr2,bbt2,bbz2)

!         implicit none
!         integer , parameter :: wp = kind(1.d0) ! working precision
!         integer , intent(in) :: lx
!         integer, intent(in)   :: m_mp, n_mp
         real(wp) , intent(in) :: amp_ar,ph_mpr
         integer(wp) , intent(in) :: n_mpr, m_mpr
         real(wp) , intent(in) :: mp_percr
!         complex(wp) , intent(out) , dimension(0:lx) :: bbr
!         complex(wp) , intent(out) , dimension(-1:lx) :: bbt, bbz
         complex(wp) , intent(out) , dimension(0:lx) :: bbr2
         complex(wp) , intent(out) , dimension(-1:lx) :: bbt2, bbz2
!         complex(wp) :: i
         real(wp)  :: dumr, dumt, dumz
         real(wp)  :: dumr2, dumt2, dumz2
!         real(wp)  :: sign_r, sign_t
!         real(wp)  :: rr
!         real(wp)  :: pi, pitol
         real(wp) , dimension(0:lx) :: axx
         real(wp) , dimension(-1:lx) :: axx2
!         integer :: ix
!         data   pi/3.14159265358979323846d0/
!         data   i/(0.d0,1.d0)/
!         data   rr/4.d0/
!         data   pitol/1.e-4/

!         write(*,*) 'prova RR module',rr,wp,lx,m_mp,n_mp,i,pi,pitol
!#! Marco, define meshes
         dx = real(1.d0/lx,wp)
         do ix=0, lx
          axx(ix) = real(ix,wp) * dx
         enddo
         do ix=-1, lx
          axx2(ix) = real((ix + 1),wp) * dx - dx/2.d0
         enddo

         sign_r = signr(ph_mpr)
         sign_t = signt(ph_mpr)
!         write(*,*) 'DIAG_apply_before',ph_mpr,dtan(ph_mpr)
!         do ix=0, lx
!           dumr = dsqrt((amp_ar * 0.5d0 
!     .           * (BESSI(m_mpr+1, n_mpr/RR*axx(IX)) 
!     .           +  BESSI(m_mpr-1,n_mpr/RR*axx(IX))))**2.d0)
!
!          if ( (abs(ph_mpr) .lt. pitol) .or. 
!     .         (abs(ph_mpr-pi) .lt. pitol) .or.
!     .         (abs(ph_mpr-2.d0*pi) .lt. pitol) ) then
!             bbr(ix) =  cmplx(0.d0 * sign_r * dumr
!     .               , (- 1.d0) * sign_r * dumr)
!          else
!             bbr(ix) = cmplx(sign_r * dumr 
!     .               , (- 1.d0) / dtan(ph_mpr) * sign_r * dumr)
!     .               * ( dsqrt( dtan(ph_mpr)**2.d0   
!     .               / (1.d0 + dtan(ph_mpr)**2.d0) ))
!          endif
!         enddo
!         write(*,*) 'DIAG_apply1',sign(1,n_mpr)
!         write(*,*) 'DIAG_apply2',sign_r * dumr
!     .              ,(- 1.d0) / dtan(ph_mpr) * sign_r * dumr

         do ix=0, lx
!          dumt = dsqrt((1.d0/axx2(IX) * amp_ar * m_mpr
!     .           * RR / n_mpr * BESSI(m_mpr, n_mpr/RR*axx2(IX)))**2.d0 )
!          dumz = dsqrt((amp_ar 
!     .           * BESSI(m_mpr, n_mpr/RR*axx2(IX)))**2.d0 )
!          if ((abs(ph_mpr-pi/2.d0) .lt. pitol) .or.
!     .        (abs(ph_mpr-3.d0*pi/2.d0) .lt. pitol)    ) then
!           bbt(ix) = cmplx( sign(1,n_mpr) * 0.d0 * sign_t * dumt
!     .            ,  sign(1,n_mpr) * sign_t * dumt )
!           bbz(ix) = cmplx(  sign_t * dumz * 0.d0
!     .            ,  sign_t * dumz )
!          else
!           bbt(ix) = cmplx( sign(1,n_mpr) * sign_t * dumt
!     .            , sign(1,n_mpr) * dtan(ph_mpr) * sign_t * dumt)
!     .            / dsqrt(1.d0 + dtan(ph_mpr)**2.d0) 
!           bbz(ix) = cmplx(  sign_t * dumz
!     .            ,  dtan(ph_mpr) * sign_t * dumz)
!     .            / dsqrt((1.d0 + dtan(ph_mpr)**2.d0)) 
!          endif
!!           bbt2(ix) = dumt * dcmplx(cos(ph_mpr+pi/2.d0)
!!     .                             ,sin(ph_mpr+pi/2.d0))
!!           bbz2(ix) = dumz * dcmplx(cos(ph_mpr+pi/2.d0)
!!     .                             ,sin(ph_mpr+pi/2.d0))
           dumt2 = 1.d0/axx2(IX) * amp_ar * m_mpr
     .           * RR / n_mpr * BESSI(m_mpr, n_mpr/RR*axx2(IX))
           dumz2 = amp_ar 
     .           * BESSI(m_mpr, n_mpr/RR*axx2(IX))
           bbt2(ix) = dumt2 * dcmplx(sin(ph_mpr)
     .                             ,-cos(ph_mpr))
           bbz2(ix) = dumz2 * dcmplx(sin(ph_mpr)
     .                             ,-cos(ph_mpr))
         enddo
!#! Marco, regularity condition
         if (m_mpr .eq. 0) then 
!          bbt(-1) = -bbt(0)
!          bbz(-1) = +bbz(0)
          bbt2(-1) = -bbt2(0)
          bbz2(-1) = +bbz2(0)
         else if (m_mpr .eq. 1) then
!          bbt(-1) = +bbt(0)
!          bbz(-1) = -bbz(0)
          bbt2(-1) = +bbt2(0)
          bbz2(-1) = -bbz2(0)
         else
!          bbt(-1) = -bbt(0)
!          bbz(-1) = -bbz(0)
          bbt2(-1) = -bbt2(0)
          bbz2(-1) = -bbz2(0)
         endif

!#! Marco, inverto i due campi perché ho un segno meno rispetto a
!PIXIE3D
!!#! Marco, 09-05-2016: provo a togliere il segno meno (rimesso)
!!#! Marco, 13-07-2016: provo a togliere il segno meno
!         bbt = - bbt
!         bbz = - bbz

!#! here I should insert a divergence cleaning part, computing br from
!bt,bz, mimando BR AUS BT +BZ
!!#! Marcp, ricorda che il br deve avere una fase pi/2 maggiore di
!quella di bt,bz
!         write(*,*) 'after'
!         bbr(0) = cmplx(0.d0,0.d0)
         bbr2(0) = cmplx(0.d0,0.d0)
         do ix=1, lx 
!          bbr(ix) = 1.d0 / axx(ix) * (
!     .              i * dx *(-m_mpr * bbt(ix-1) 
!     .              - n_mpr / rr * axx2(ix-1) * bbz(ix-1) ) 
!     .              + bbr(ix-1) * axx(ix-1) )
          bbr2(ix) = 1.d0 / axx(ix) * (
     .              i * dx *(-m_mpr * bbt2(ix-1) 
     .              - n_mpr / rr * axx2(ix-1) * bbz2(ix-1) ) 
     .              + bbr2(ix-1) * axx(ix-1) )
         enddo
         if (m_mpr .eq. 1) then
!          bbr(0) = (4.d0*bbr(1) - bbr(2)) / 3.d0
!          bbt(-1) = - bbt(0) +2.d0*i*bbr(0)
          bbr2(0) = (4.d0*bbr2(1) - bbr2(2)) / 3.d0
          bbt2(-1) = - bbt2(0) +2.d0*i*bbr2(0)
         endif
!         do ix=0,lx 
!          write(*,*) 'DIAG ix=',ix
!          write(*,*) bbr(ix),bbr2(ix)
!          write(*,*) bbt(ix),bbt2(ix)
!          write(*,*) bbz(ix),bbz2(ix)
!         enddo

!!#! Marco, check wheter the amplitude at edge is correct
!#! Marco, not necessary
!!          write(*,*) 'DIAGcrt_fct',mp_percr
!          crt_fct = mp_percr / abs(bbr(lx))
!!          write(*,*) 'DIAGcrt_fct',crt_fct
!          bbr = bbr * crt_fct
!          bbt = bbt * crt_fct
!          bbz = bbz * crt_fct
!
!         write(*,*) 'DIAG_apply2',bbr
        end subroutine apply_irrot_bc

      function signr(phmpr)
!#! Marco, building signs of the real and imaginary parts of magnetic
!field components, depending on the value of ph_mp
!#!Marco, cambiato il 14/07/2016
       implicit none
       real (wp) , intent(in) :: phmpr
       real (wp) :: signr
        if (phmpr .ge. 0.d0 .and. phmpr .lt. pi/2.d0) then
         signr = - 1. 
        else if (phmpr .ge. pi/2.d0 .and. phmpr .lt. pi) then
         signr = - 1. 
        else if (phmpr .ge. pi .and. phmpr .lt. 3.d0*pi/2.d0) then
         signr = + 1. 
        else if (phmpr .ge. 3.d0*pi/2.d0 .and. phmpr .lt. 2.d0*pi) then
         signr = + 1. 
        endif
      end function signr
       
      function signt(phmpr)
!#! Marco, building signs of the real and imaginary parts of magnetic
!field components, depending on the value of phmpr
       implicit none
       real (wp) , intent(in) :: phmpr
       real (wp) :: signt
        if (phmpr .ge. 0.d0 .and. phmpr .lt. pi/2.d0) then
         signt = + 1. !* -1.
        else if (phmpr .ge. pi/2.d0 .and. phmpr .lt. pi) then
         signt = - 1. ! * -1.
        else if (phmpr .ge. pi .and. phmpr .lt. 3.d0*pi/2.d0) then
         signt = - 1.  !* -1.
        else if (phmpr .ge. 3.d0*pi/2.d0 .and. phmpr .lt. 2.d0*pi) then
         signt = + 1.  !* -1.
        endif
      end function signt

      function get_phase(field)
!#! Marco, get the correct phase from atan
       implicit none
       complex (8) , intent(in) :: field
       real (8) :: get_phase, pi, phase
       DATA   pi/3.14159265358979323846d0/

       phase = atan( aimag(field) / real(field) )

       if (aimag(field) .ge. 0.d0) then
        if (real(field) .ge. 0d0) then
         get_phase = phase         
        else
!#! Marco, phase in questo caso è negativa
         get_phase = phase + pi
        endif
       else 
        if (real(field) .ge. 0d0) then
!#! Marco, phase in questo caso è negativa
         get_phase = 2.d0*pi + phase
        else 
         get_phase = phase + pi
        endif
       endif

      end function get_phase

      subroutine apply_mp_init
!#! Marco, subroutine that applies the magnetic field boundary
!conditions at the beginning of each part of the simulation

      !local variables
      logical :: log1, log2

! Daniele e Marco, 02/10/2014
! inizializziamo br con profilo di simil vuoto
! e bt accordingly per avere div(B)=0
! SpeCyl se lo mangia molto bene
! In futuro bisognerebbe sistemare meglio il codice
! e implementare il caso Bessel con anche rot(B)=0
! migliore per il caso stellarator
! e in vista di variazione dinamica delle BCs
!#! Marco, definition of some varibles needed also on the case

!#! inizioparteMP_Marco

      do idx=0, size(n_mp,1)-1 
       if (m_mp(idx) .eq. 0 .and. n_mp(idx) .eq. 0
     .     .and. mp_perc(idx) .eq. 0.d0) then
        write(*,*) 'non ci sono MP da imporre, esco e vado'
        exit
       endif
       write(*,*) 'mp0',BESSI(m_mp(idx)+1, n_mp(idx)/RR)
     .              +  BESSI(m_mp(idx)-1, n_mp(idx)/RR)
!       write(*,*) m_mp(idx),n_mp(idx),janz(m_mp(idx),n_mp(idx)),RR
       j_mp(idx) = janz(m_mp(idx),n_mp(idx))
       j_mpold(idx) = janz(m_mpold(idx),n_mpold(idx))
       write(*,'(a,4i3,f5.2)') 'check',m_mp(idx)
     .         ,n_mp(idx),j_mp(idx),j_mpold(idx),RR
       if (mp_perc(idx) .ne. 0.d0) then 
        amp_a(idx) = 2.d0 * mp_perc(idx) /
     .               (BESSI(m_mp(idx)+1, n_mp(idx)/RR)
     .                +  BESSI(m_mp(idx)-1, n_mp(idx)/RR))
!         write(*,*) 'ampa',amp_a(idx)
       endif
       write(*,*) 'inizioMP2' ,amp_a(idx)
       call  flush(6)
       if (ph_mp(idx) .lt. 0.d0 .and. ph_mp(idx) .gt. 2.d0*pi) then
       write(*,*) 'Use a phase between zero and two pi', ph_mp(idx)
       stop
       endif
       !#! Marco, 03/12/2014 implemento anche le MP irrotazionali
!       write(*,*) 'DEBUG-1',j_mp,j_mpold
       call flush(6)
       if (.not. irrot) then !(@1)
        IF (.NOT.IBAND) THEN !(@2)
!#! Marco, 2018-02-11 inserire qui le MP rotazionali per i casi tokamak
!con profilo r^(m+tkp)
         if (tok_rot) then !#!(@3)
          if (tau_growth .lt. 1.d0) then !#! @4!in this case I apply the MP at
!#! the beginning
           write(*,*) 'perturbing tokamak simulation with tkp=',tkp
           write(*,*) 'mp_perc=',mp_perc(idx)
!          IF (amp_a(idx) .eq. 0.d0) cycle
           DO IX=0,LX
            IF (N.EQ.0) THEN
             BR(IX,j_mp(idx))=mp_perc(idx)*I*X(IX)**(m_mp(idx)-1)
            ELSE
! perturbare solo br e bt come in PIXIE3D toroidale, br=r^(m+tkp)
             BR(IX,j_mp(idx))=mp_perc(idx)*I*X(IX)**(m_mp(idx)+tkp)
            END IF
            BR(IX,-j_mp(idx))=CONJG(BR(IX,j_mp(idx)))

!                   write(6,*)'r,br ',X(IX),AIMAG(BR(IX,J))
!                   call flush(6)
           END DO
           DO IX=0,LX
            IF (N.EQ.0) THEN
             BT(IX,j_mp(idx))=-mp_perc(idx)*X2(IX)**(m_mp(idx)-1)
             BZ(IX,j_mp(idx))=0.
            ELSE
! modo di perturbare solo br e bt come in PIXIE3D toroidale, br=r^(m+tkp)
             BT(IX,j_mp(idx))=-mp_perc(idx)*(m_mp(idx)+tkp+1)
     .                       /m_mp(idx)*X2(IX)**(m_mp(idx)+tkp)
             BZ(IX,j_mp(idx))=0.
            END IF
            BT(IX,-j_mp(idx))=CONJG(BT(IX,j_mp(idx)))
            BZ(IX,-j_mp(idx))=CONJG(BZ(IX,j_mp(idx)))
           END DO
           BT(-1,j_mp(idx))=-BT(1,j_mp(idx))
           BT(-1,-j_mp(idx))=CONJG(BT(-1,j_mp(idx)))
           BZ(-1,j_mp(idx))=-BZ(1,j_mp(idx))
           BZ(-1,-j_mp(idx))=CONJG(BZ(-1,j_mp(idx)))
          endif !#!@4
         else !#!@3 tok_rot

!         write(*,*) 'implement rotational MPs!'
           J=1
           M=1
!           AMPA=-1.5E-3
           AMPB=3.
           DO IX=0,LX
            BR(IX,J)=mp_perc(idx)*I*X(IX)**AMPB
            BR(IX,-J)=CONJG(BR(IX,J))
           END DO
           DO IX=-1,LX
            BT(IX,J)=-mp_perc(idx)/M*(AMPB+1)*X2(IX)**AMPB
            BT(IX,-J)=CONJG(BT(IX,J))
            BZ(IX,J)=0.
            BZ(IX,-J)=CONJG(BZ(IX,J))
           END DO
          endif!(@3) fine caso rotazionale
         else !#! (@2) inizio caso in cui si ricomincia la simulazione
          !#! non dovrebbe servire nulla qui 
         endif !#! (@2) fine caso in cui si ricomincia la simulazione
        else !(@1)
       !#! irrotational MP
       !#! Marco, 19/01/15: if simulation is in beginning mode APPLY the
       !irrotational MP of the selected helicity
        if (.not. IBAND) then !(@2)
!        write(*,*) '#*#*#*', mp_perc(idx)
!        write(*,*) '#@', BESSI(m_mp(idx)+1, n_mp(idx)/RR)
        write(*,*) 'applying irrotational MP on modes('
     .               ,m_mp(idx),',',n_mp(idx),') ',j_mp(idx)
        write(*,*) 'applying irrotational MP with amp and phase'
        write(*,*) mp_perc(idx), ph_mp(idx), 'normalized amplitude'
        write(*,*) 'applying irrotational MP rotating with freq'
        write(*,*) br_pert_freq
        if (j_mpold(idx) /= 0) then
         write(*,*) 'coming from MP on modes('
     .               ,m_mpold(idx),',',n_mpold(idx),') ',j_mpold(idx)
        endif
        call flush(6)
         if (mp_perc(idx) .ne. 0.d0) then !(@3)
          amp_a(idx) = 2.d0 * mp_perc(idx) 
     .                 / (BESSI(m_mp(idx)+1, n_mp(idx)/RR)
     .                 +  BESSI(m_mp(idx)-1, n_mp(idx)/RR))
          write(*,*) 'INIZIO LA SIMULAZIONE E APPLICO MP IRROT'
          write(*,*) 'amp_a', amp_a(idx)
          call apply_irrot_bc(mp_perc(idx),amp_a(idx),m_mp(idx)
     .             ,n_mp(idx),ph_mp(idx),bbr,bbt,bbz)
          do ix=0, lx 
           br(ix,j_mp(idx)) = bbr(ix)
           br(ix,-j_mp(idx)) = conjg(br(ix,j_mp(idx)))
          enddo
          do ix=-1, lx 
           bt(ix,j_mp(idx)) = bbt(ix)
           bz(ix,j_mp(idx)) = bbz(ix)
           bt(ix,-j_mp(idx)) = conjg(bt(ix,j_mp(idx)))
           bz(ix,-j_mp(idx)) = conjg(bz(ix,j_mp(idx)))
          enddo
!          write(*,*) 'check amplitude', abs(bbr(lx)),abs(br(lx,73))
       !#! Marco, end of mp_perc ne 0
         endif !(@3)
       !#! Marco, inizializzo la variabile phmp, utile per rotating MP
         phmp(idx) = ph_mp(idx)
       !#! Marco, if br_pert_freq /= 0 leggo l'ultima frequenza
       !#! This is not needed, if the simulation is starting there is no need
       !to read the old frequency
       !          if (br_pert_freq /= 0.d0) then
       !           phmp(idx) = atan( aimag(br(lx,j_mp))/real(br(lx,j_mp)) ) + pi/2.d0
       !           write(*,*) 'Rotating MP, the old frequency is', phmp(idx)
       !          endif
       !#! Marco, 19/01/15: if simulation is in restart mode SUBTRACT the old

       !irrotational MP and ADD the new ones
        else !(@2) se sto continuando la simulazione
         if (mp_perc(idx) .ne. 0.d0) then !(@3)
         write(*,*) 'CONTINUO LA SIMULAZIONE E APPLICO MP IRROT'
         call flush(6)
       !#! Marco , if br_pert_freq /= 0 leggo l'ultima frequenza (dall'input
       !file)
         if (br_pert_freq /= 0.d0) then !(@4)
          if (ampold(idx) .eq. 0.d0) then !(@5)
           write(*,*) 'Rotating MP: starting from nonMP case'
     . ,                amp_a(idx),ph_mp(idx)
       !#! Marco, add new irrotational perturbation 
           call apply_irrot_bc(mp_perc(idx),amp_a(idx),m_mp(idx)
     .          ,n_mp(idx),ph_mp(idx),bbr,bbt,bbz)
           do ix=0, lx 
            br(ix,j_mp(idx)) = br(ix,j_mp(idx)) + bbr(ix)
            br(ix,-j_mp(idx)) = conjg(br(ix,j_mp(idx)))
           enddo
           do ix=-1, lx 
            bt(ix,j_mp(idx)) = bt(ix,j_mp(idx)) + bbt(ix)
            bz(ix,j_mp(idx)) = bz(ix,j_mp(idx)) + bbz(ix)
            bt(ix,-j_mp(idx)) = conjg(bt(ix,j_mp(idx)))
            bz(ix,-j_mp(idx)) = conjg(bz(ix,j_mp(idx)))
           enddo
       !#! Marco , inizializzo la variabile phmp, utile per rotating MP
           phmp(idx) = ph_mp(idx)
          else !(@5)
           write(*,*) 'Rotating MP: starting from a MP case. I need to
     . get the old phase and amplitude'
           
!           phmp(idx) = atan( aimag(br(lx,j_mp(idx)))
!     .                 /real(br(lx,j_mp(idx))) ) 
           phmp(idx) = get_phase(br(lx,j_mp(idx)))
!       !#! Marco (spiegazione nuova, del 14/07/2016) metto -pi/2.d0 perché la fase è
!       quella delle componenti th,z, che è -pi/2. rispetto a quella del
!       br
     .                - pi/2.d0
!        write(*,*) 'DIAG rst get_phase bz',phmp(idx),'
!     .  br',phmp(idx)+pi/2.d0
!        write(*,*) 'DIAG restart get_br0',br(lx,j_mp(idx))
!       !#! Marco (spiegazione vecchia), metto -pi/2.d0 perché nella ruotine apply_irrot_bc impongo
!       !d'autorità un altro - alle componenti theta e zeta e quindi la fase del
!       !theta, che sarebbe più pi/2 rispetto a quella della componente r
!       !diventa meno pi/2
           if (phmp(idx) .lt. 0.d0) then 
            phmp(idx) = phmp(idx) + 2.d0*pi
           endif
        !#! Marco, a questo punto devo sottrarre le vecchie MP e
        !aggiungere quelle nuove
           amp_old = 2.d0 * ampold(idx) 
     .               / (BESSI(m_mpold(idx)+1, n_mpold(idx)/RR)
     .               +  BESSI(m_mpold(idx)-1, n_mpold(idx)/RR))
           call apply_irrot_bc(ampold(idx),amp_old,m_mpold(idx),
     .                           n_mpold(idx),phmp(idx),bbr,bbt,bbz)
!           write(*,*) 'DIAG restart get_br1',bbr(lx),bbz(lx)
           do ix=0, lx 
            br(ix,j_mpold(idx)) = br(ix,j_mpold(idx)) - bbr(ix)
            br(ix,-j_mpold(idx)) = conjg(br(ix,j_mpold(idx)))
           enddo
           do ix=-1, lx 
            bt(ix,j_mpold(idx)) = bt(ix,j_mpold(idx)) - bbt(ix)
            bz(ix,j_mpold(idx)) = bz(ix,j_mpold(idx)) - bbz(ix)
            bt(ix,-j_mpold(idx)) = conjg(bt(ix,j_mpold(idx)))
            bz(ix,-j_mpold(idx)) = conjg(bz(ix,j_mpold(idx)))
           enddo

       !#! Marco, add new irrotational perturbation 
           call apply_irrot_bc(mp_perc(idx),amp_a(idx),m_mp(idx)
     .         ,n_mp(idx),phmp(idx),bbr,bbt,bbz)
!           write(*,*) 'add new MP'
!           write(*,*) 'DIAG restart get_br2',bbr(lx),bbz(lx)
           do ix=0, lx 
            br(ix,j_mp(idx)) = br(ix,j_mp(idx)) + bbr(ix)
            br(ix,-j_mp(idx)) = conjg(br(ix,j_mp(idx)))
           enddo
           do ix=-1, lx 
            bt(ix,j_mp(idx)) = bt(ix,j_mp(idx)) + bbt(ix)
            bz(ix,j_mp(idx)) = bz(ix,j_mp(idx)) + bbz(ix)
            bt(ix,-j_mp(idx)) = conjg(bt(ix,j_mp(idx)))
            bz(ix,-j_mp(idx)) = conjg(bz(ix,j_mp(idx)))
           enddo
!           phmp(idx) = atan2(aimag(br(lx,j_mp(idx))),
!     .                      real(br(lx,j_mp)))
           phmp(idx) = get_phase(br(lx,j_mp(idx)))
!           write(*,*) 'after sub - add', real(br(lx,j_mp(idx)))
!     .                     , aimag(br(lx,j_mp(idx)))
!     .                   , phmp(idx)
          endif  !(@5)
         else !(@4)!#! Marco, parte in cui continuo da una simulazione precedente ma br_pert_freq=0
       !#! ampold(idx) può essere calcolato già in base ai dati che
       !abbiamo letto dal file di input
           ampold(idx) = abs(br(lx,j_mp(idx)))
           if (ampold(idx) .eq. 0.) then 
            phold = 0.d0
            amp_old = 0.d0
           else
            phold = get_phase(br(lx,j_mp(idx))) - pi/2.d0
            amp_old = 2.d0 * ampold(idx) 
     .               / (BESSI(m_mpold(idx)+1, n_mpold(idx)/RR)
     .               +  BESSI(m_mpold(idx)-1, n_mpold(idx)/RR))
           endif
           write(*,*) 'old phase',phold,amp_old,ampold,mp_perc(idx)
           log1 = amp_old .eq. 0.d0
!           log2 = abs(amp_old - mp_perc(idx)) .lt. pitol !#! cambio
!           questa cosa, perché amp_old è già convertita in bessel
!           mentre me_perc è quella che do io in input. Uso piuttosto
!           ampold
           log2 = abs(ampold(idx) - mp_perc(idx)) .lt. pitol
           write(*,*) 'logical',log1,log2
           if (log1 .or. log2 ) then !(@5)
            write(*,*) 'nothing to subtract / nothing to do'
           else !(@5)
            call apply_irrot_bc(ampold(idx),amp_old,m_mpold(idx),
     .                           n_mpold(idx),phold,bbr,bbt,bbz)
            write(*,*) 'subtract old MP', j_mpold(idx)
            do ix=0, lx 
             br(ix,j_mpold(idx)) = br(ix,j_mpold(idx)) - bbr(ix)
             br(ix,-j_mpold(idx)) = conjg(br(ix,j_mpold(idx)))
            enddo
            do ix=-1, lx 
             bt(ix,j_mpold(idx)) = bt(ix,j_mpold(idx)) - bbt(ix)
             bz(ix,j_mpold(idx)) = bz(ix,j_mpold(idx)) - bbz(ix)
             bt(ix,-j_mpold(idx)) = conjg(bt(ix,j_mpold(idx)))
             bz(ix,-j_mpold(idx)) = conjg(bz(ix,j_mpold(idx)))
            enddo
       
       !#! Marco, add new irrotational perturbation 
            ph_mp(idx) = get_phase(br(lx,j_mp(idx))) - pi/2.d0
            call apply_irrot_bc(mp_perc(idx),amp_a(idx),m_mp(idx)
     .          ,n_mp(idx),ph_mp(idx),bbr,bbt,bbz)
            write(*,*) 'add new MP, br_pert_freq=0'
            do ix=0, lx 
             br(ix,j_mp(idx)) = br(ix,j_mp(idx)) + bbr(ix)
             br(ix,-j_mp(idx)) = conjg(br(ix,j_mp(idx)))
            enddo
            do ix=-1, lx 
             bt(ix,j_mp(idx)) = bt(ix,j_mp(idx)) + bbt(ix)
             bz(ix,j_mp(idx)) = bz(ix,j_mp(idx)) + bbz(ix)
             bt(ix,-j_mp(idx)) = conjg(bt(ix,j_mp(idx)))
             bz(ix,-j_mp(idx)) = conjg(bz(ix,j_mp(idx)))
            enddo
           endif !(@5)
       !#! Marco , inizializzo la variabile phmp, utile per rotating MP
           phmp(idx) = ph_mp(idx)
           write(*,*) 'is defined',phmp(idx)
       !#! Marco, end of the br_pert_freq
         endif !(@4)
       !#! Marco, end of the mp_perc ne 0 clause
        endif !(@3)
       !#! Marco, end of the IBAND clause
       endif !(@2)
       !#! Marco, end of the irrot =T or F clause
       endif !(@1)
       !#! Marco, end of the cycle on n_mp
       enddo
       call flush(6)

       write(*,*)  'applied MP',j_mp(0), br(lx,j_mp(0))
     .           ,bt(lx,j_mp(0)), bz(lx,j_mp(0))

!#! fineparte

      end subroutine apply_mp_init

      subroutine apply_mp
!#! Marco, subroutine that applies the magnetic field boundary
!conditions 
!#! Marco, attenzione per MP irrotazionali
!#! Marco, escludo anche il caso di MP rotanti, l'applicazione delle BC
!viene fatta, in quel caso, più avanti nel codice 
         if (irrot .eqv. .true. .and. mp_perc(0) .ne. 0.d0) then !#!@1
          if (br_pert_freq .eq. 0.d0) then !#! @2
!#! remind, defined in apply_mp_init
!       amp_a(idx) = 2.d0 * mp_perc(idx) 
!     .              / (BESSI(m_mp(idx)+1, n_mp(idx)/RR)
!     .                 +  BESSI(m_mp(idx)-1, n_mp(idx)/RR))
           call flush(6)
!#! inizioparte3
           if (any(j_mp .eq. j)) then !#!@3
            do ik=0, qmp
             temp_mp(ik) = abs(j_mp(ik)-j)
            end do
            loc = minloc(temp_mp,1) - 1
            mloc = m_mp(loc)
            nloc = n_mp(loc)
            sign_r = signr(phmp(loc))
!#! Marco, 2 novembre 2018. Inserisco delle MP che aumentano nel tempo,
!usando le keywords mp_exp,mp_tanh
!#! Marco, qui imposto solo il valore della MP al bordo
             if (mp_exp) then
              UU(1,lx) = mp_cmplx(loc)
     .                 * (heaviside(time_passage-tot_it*dt) 
     .                 * (1.d0-exp(-((tot_it-it_beginning)*dt)
     .                 /tau_growth)) 
     .                 + heaviside((tot_it+1)*dt-time_passage)
     .                 *exp(-(tot_it*dt-time_passage)/tau_decrease) )

!              if (j .eq. 74 .and. (tot_it mod 100 .eq. 0)) then 
!               write(*,*) 'DIAG mp_exp',UU(1,lx),time_passage
!     .                    ,tau_decrease
!     .                    ,heaviside((tot_it+1)*dt-time_passage)
!     .                    ,(tot_it)*dt
!              endif
             endif
             if (mp_tanh) then
              UU(1,lx) = mp_cmplx(loc)/2.d0
     .            * ( heaviside(time_passage-tot_it*dt) 
     .            * (tanh((tot_it*dt-time_medium)/tau_growth) + 1.d0) 
     .            + heaviside((tot_it+1)*dt-time_passage)
     .            * (tanh(-(tot_it*dt - time_medium - time_passage)
     .            /tau_decrease)+ 1.d0 ) )
             endif

!            UU(1,lx) = cmplx(real(mp_cmplx(0)),aimag(mp_cmplx(0)))

!            dumr = sqrt(( amp_a(loc) * 0.5d0 
!     .           * (BESSI(mloc+1, nloc/RR)
!     .           +  BESSI(mloc-1, nloc/RR) ) )**2.d0)
!
!             if ( (abs(phmp(loc)) .lt. pitol) .or. 
!     .            (abs(phmp(loc)-pi) .lt. pitol) .or.
!     .            (abs(phmp(loc)-2.d0*pi) .lt. pitol) ) then
!                UU(1,lx) =  cmplx(0.d0 * sign_r * dumr
!     .                  , (- 1.d0) * sign_r * dumr)
!             else
!                UU(1,lx) = cmplx(sign_r * dumr 
!     .            ,  (- 1.d0) / tan(phmp(loc)) * sign_r * dumr)
!     .            * ( sqrt( tan(phmp(loc))**2.d0   
!     .            / (1.d0 + tan(phmp(loc))**2.d0) ))
!             endif
!#! fineparte3
!#! Marco, 13 luglio 2017. Introduco l'input mp_cmplx, per non dover
!fare la conversione modulo-fase --> numero complesso che è irta di
!ambiguità
            UU(1,lx) = mp_cmplx(loc)
!            write(*,*) 'inside_mp, check UU',UU(1,lx)
!          write(*,*) 'inside_mp, check UU / second',UU(1,lx)
           endif !#!@3
          else !#!@2 (br_pert_freq /= 0.d0)
!#! attenzione, solo se it = 1
           call flush(6)
!!#! inizio (br_pert_freq /= 0)
           if (any(j_mp .eq. j)) then
!#! Marco, attenzione. Quando impongo la condizione al contorno devo
!ricordare che phmp(loc) è la fase di bz,bt!
           write(*,*) 'qui solo se le MP ruotano'
           do ik=0, qmp
            temp_mp(ik) = abs(j_mp(ik)-j)
           end do
           loc = minloc(temp_mp,1) - 1
           mloc = m_mp(loc)
           nloc = n_mp(loc)
           sign_r = signr(phmp(loc))
!           write(*,*) 'apply_bc_final',loc,phmp(loc)
           dumr = sqrt(( amp_a(loc) * 0.5d0 
     .              * (BESSI(mloc+1, nloc/RR)
     .              +  BESSI(mloc-1, nloc/RR) ) )**2.d0)

           UU(1,lx) = dumr * cmplx(cos(phmp(loc)),sin(phmp(loc)))
!           if (j .eq. 1) then 
!            write(*,*) 'DIAGapp',phmp(loc),UU(1,lx)
!           endif
!            if ( (abs(phmp(loc) ) .lt. pitol) .or. 
!     .           (abs(phmp(loc) - pi) .lt. pitol) .or.
!     .           (abs(phmp(loc)- 2.d0*pi) .lt. pitol) ) then
!               UU(1,lx) =  cmplx(0.d0 * sign_r * dumr
!     .                 , (- 1.d0) * sign_r * dumr)
!            else
!               UU(1,lx) = cmplx(sign_r * dumr 
!     .           ,  (- 1.d0) / tan(phmp(loc)) * sign_r * dumr)
!     .           * ( sqrt( tan(phmp(loc))**2.d0   
!     .           / (1.d0 + tan(phmp(loc))**2.d0) ))
!            endif
!#! fineparte3
           endif
          
          endif!#!@2 (end of br_pert_freq /= 0.d0)
!         endif!#!@1
          else !#!@1 (here case for rotational MPs)
!           write(*,*) 'qui solo se MP rotational'
           if (mp_perc(0) .ne. 0) then
            if (tok_rot) then
!#! logical check on mpexp,mptanh
             if (mptanh) then
              mpexp = .false.
             endif
             if (br_pert_freq .eq. 0.d0) then !#! @2
              if (any(j_mp .eq. j)) then !#!@3
!               write(*,*) 'sto imponendo MP rotazionali sul caso
!     .         tokamak',j, mp_cmplx,br(lx,j),mp_perc(0),tot_it
!               UU(1,lx) = mp_cmplx(0)
!               UU(1,lx) = br(lx,j)

               if (mpexp) then
                UU(1,lx) = mp_cmplx(0)
     .                   * (heaviside(time_passage-tot_it*dt) 
     .                   * (1.d0-exp(-(tot_it*dt)/tau_growth)) 
     .                   + heaviside((tot_it+1)*dt-time_passage)
     .                   *exp(-(tot_it*dt-time_passage)/tau_decrease) )
               endif
               if (mptanh) then
                UU(1,lx) = mp_cmplx(0)/2.d0
     .              * ( heaviside(time_passage-tot_it*dt) 
     .              * (tanh((tot_it*dt-time_medium)/tau_growth) + 1.d0) 
     .              + heaviside((tot_it+1)*dt-time_passage)
     .              * (tanh(-(tot_it*dt - time_medium - time_passage)
     .              /tau_decrease)+ 1.d0 ) )
               endif
!               write(*,*) UU(1,lx)
              endif
             else !#! MP rotazionali e rotanti sul tokamak con
              if (any(j_mp .eq. j)) then !#!@3
               do ik=0, qmp
                temp_mp(ik) = abs(j_mp(ik)-j)
               end do
               loc = minloc(temp_mp,1) - 1
!               UU(1,lx) = abs(mp_cmplx(0)) 
!     .                  * cmplx(cos(phmp(loc)),sin(phmp(loc)))
!#! Marco, boundary condition modulated in amplitude and phase
               if (mpexp) then
               UU(1,lx) = mp_perc(loc) 
     .                  * cmplx(cos(phmp(loc)),sin(phmp(loc)))
     .                  * (heaviside(time_passage-tot_it*dt) 
     .                  * (1.d0-exp(-(tot_it*dt)/tau_growth)) 
     .                  + heaviside((tot_it+1)*dt-time_passage)
     .                  *exp(-(tot_it*dt-time_passage)/tau_decrease) )
               endif
               if (mptanh) then
                UU(1,lx) = mp_cmplx(0)/2.d0
     .               * cmplx(cos(phmp(loc)),sin(phmp(loc)))
     .               * ( heaviside(time_passage-tot_it*dt) 
     .               * (tanh((tot_it*dt-time_medium)/tau_growth) + 1.d0) 
     .               + heaviside((tot_it+1)*dt-time_passage)
     .               * (tanh(-(tot_it*dt - time_medium - time_passage)
     .                 /tau_decrease)) + 1.d0 )
               endif
              endif
             endif
            endif
           endif
          endif

      end subroutine apply_mp

!************************************************************************
!*                                                                      *
!*    Program to calculate the first kind modified Bessel function of   *
!*    integer order N, for any REAL X, using the function BESSI(N,X).   *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*                                                                      *
!*    SAMPLE RUN:                                                       *
!*                                                                      *
!*    (Calculate Bessel function for N=2, X=0.75).                      *
!*                                                                      *
!*    Bessel function of order  2 for X =  0.7500:                      *
!*                                                                      *
!*         Y =  0.73666878E-01                                          *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!*                                                                      *
!*                               F90 Release 1.1 By J-P Moreau, Paris.  *
!*                                        (www.jpmoreau.fr)             *
!*                                                                      *
!*   Version 1.1: corected value of P4 in BESSIO (P4=1.2067492 and not  *
!*                1.2067429) Aug. 2011.                                 *
!************************************************************************
c$$$PROGRAM TBESSI
c$$$
c$$$  REAL*8  BESSI, X, Y
c$$$  INTEGER N
c$$$
c$$$  N=2
c$$$  X=0.75D0
c$$$
c$$$  Y = BESSI(N,X)
c$$$
c$$$  write(*,10) N, X
c$$$  write(*,20) Y
c$$$
c$$$  stop
c$$$
c$$$ 10    format (/' Bessel Function of order ',I2,' for X=',F8.4,':')
c$$$ 20     format(/'      Y = ',E15.8/)
c$$$
c$$$END

! ----------------------------------------------------------------------
      FUNCTION BESSI(NBBESS,XBBESS)
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order NBBESS, for any REAL XBBESS. We use here the classical
!     recursion formula, when XBBESS > NBBESS. For XBBESS < NBBESS, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      integer(wp) ,intent(in) :: nbbess
      REAL(wp) :: XBBESS,BESSI,TOX,BIM,BI,BIP
      IF (NBBESS.EQ.0) THEN
      BESSI = BESSI0(XBBESS)
      RETURN
      ENDIF
      IF (NBBESS.EQ.1) THEN
      BESSI = BESSI1(XBBESS)
      RETURN
      ENDIF
      IF(XBBESS.EQ.0.D0) THEN
      BESSI=0.D0
      RETURN
      ENDIF
      TOX = 2.D0/XBBESS
      BIP = 0.D0
      BI  = 1.D0
      BESSI = 0.D0
      MBBESS = 2*((NBBESS+INT(SQRT(FLOAT(IACC*NBBESS)))))
      DO 12 JBBESS = MBBESS,1,-1
      BIM = BIP+DFLOAT(JBBESS)*TOX*BI
      BIP = BI
      BI  = BIM
      IF (ABS(BI).GT.BIGNO) THEN
      BI  = BI*BIGNI
      BIP = BIP*BIGNI
      BESSI = BESSI*BIGNI
      ENDIF
      IF (JBBESS.EQ.NBBESS) BESSI = BIP
 12    CONTINUE
      BESSI = BESSI*BESSI0(XBBESS)/BI
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------
! Auxiliary Bessel functions for N=0, N=1
      FUNCTION BESSI0(XBBESS)
      REAL(8) :: XBBESS,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,
     .     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,
     .     1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     .     0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     .     0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(XBBESS).LT.3.75D0) THEN
      Y=(XBBESS/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(XBBESS)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------
      FUNCTION BESSI1(XBBESS)
      REAL(8) :: XBBESS,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,
     .     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     .     0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     .     -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     .     -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(XBBESS).LT.3.75D0) THEN
      Y=(XBBESS/3.75D0)**2
      BESSI1=XBBESS*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(XBBESS)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------

      end module mp
