      module equations_solution

      use defs
!#! Marco, subroutines and functions to solve the visco-resistive MHD
!equations in cylindrical geometry
      contains

      subroutine faraday(iilx,iiiz,br_in,bt_in,bz_in,vr_in,vt_in
     .                  ,vz_in,br_out,bt_out,bz_out,d_t,sseq)
  
       integer(wp) , intent(in) :: iilx, iiiz, sseq
       complex(wp) , intent(in) :: br_in(0:iilx,-iiiz:iiiz)
       complex(wp) , intent(in) :: bt_in(-1:iilx,-iiiz:iiiz)
       complex(wp) , intent(in) :: bz_in(-1:iilx,-iiiz:iiiz)
       complex(wp) , intent(in) :: vr_in(0:iilx,-iiiz:iiiz)
       complex(wp) , intent(in) :: vt_in(-1:iilx,-iiiz:iiiz)
       complex(wp) , intent(in) :: vz_in(-1:iilx,-iiiz:iiiz)
       complex(wp) , intent(out) :: br_out(0:iilx,-iiiz:iiiz)
       complex(wp) , intent(out) :: bt_out(-1:iilx,-iiiz:iiiz)
       complex(wp) , intent(out) :: bz_out(-1:iilx,-iiiz:iiiz)
!       complex(wp) , intent(in) :: br_in(:,:)
!       complex(wp) , intent(in) :: bt_in(:,:)
!       complex(wp) , intent(in) :: bz_in(:,:)
!       complex(wp) , intent(in) :: vr_in(:,:)
!       complex(wp) , intent(in) :: vt_in(:,:)
!       complex(wp) , intent(in) :: vz_in(:,:)
!       complex(wp) , intent(out) :: br_out(:,:)
!       complex(wp) , intent(out) :: bt_out(:,:)
!       complex(wp) , intent(out) :: bz_out(:,:)
       real(wp) , intent(in) :: d_t

! Local variables
  

! computational part

       DO IM=0,MY       
       M = MM(IM)          
       DO N=NA1(M),NE(M)
         SUM0=(0.0d0,0.0d0)
         SUM1=(0.0d0,0.0d0)
         SUM2=(0.0d0,0.0d0)
         SUM3=(0.0d0,0.0d0)
         SUM8=(0.0d0,0.0d0)
       
       nj = 0
       jmn = janz(m,n)
!       write(*,*) 'DIAGfar',vz_in(10,0),vt_in(10,1)
!       write(*,*) 'DIAGfarad-11',bz_in(:,0)
!       write(*,*) 'DIAGfarad-2',bz_in(:,1)
!       write(*,*) 'DIAGfarad-11',vz_in(:,0)
!       write(*,*) 'DIAGfarad-2',j,br_in(:,j)
!       write(*,*) 'DIAGfarad-2',-j,br_in(:,-j)
!       write(*,*) 'DIAGfarad-3',j,vt_in(:,j)
!       write(*,*) 'DIAGfarad-4',m,n,convnj(jmn),na1(m),ne(m)
       if (m.eq.0) jmn = -jmn ! serve altrimenti ho jmn negativi!!!
!$OMP PARALLEL DO DEFAULT(SHARED)
!#! Marco
!$OMP& PRIVATE(IAM,js,IX,J,J1,m,n)
!$OMP& SCHEDULE(STATIC)
!$OMP& REDUCTION(+:SUM0,SUM1,SUM2,SUM3,SUM8)
       DO js=0,convnj(jmn)-1
       IAM=OMP_GET_THREAD_NUM()
       j = convj(jmn,js)
       j1 = convj1(jmn,js)
!       if (jmn .eq. 0) then
!        write(*,'(a,4i4)') 'aaa',js,jmn,j,j1
!       endif
       DO IX=0,LX-1                                                  
       SUM3(IX) = SUM3(IX)+vt_in(IX,J)*bz_in(IX,J1)
     .                    -vz_in(IX,J)*bt_in(IX,J1)
       enddo
       DO IX=0,LX                                                    
        SUM0(IX) = SUM0(IX)+br_in(IX,J)*(vt_in(IX-1,J1)+vt_in(IX,J1))
        SUM1(IX) = SUM1(IX)+vr_in(IX,J)*(bt_in(IX-1,J1)+bt_in(IX,J1))
        SUM2(IX) = SUM2(IX)+br_in(IX,J)*(vz_in(IX-1,J1)+vz_in(IX,J1))
        SUM8(IX) = SUM8(IX)+vr_in(IX,J)*(bz_in(IX-1,J1)+bz_in(IX,J1))
!        write(*,*) 'DF-4',ix,j,j1,sum0(ix),vt_in(ix,j1),vz_in(ix,j1)
!     .                   ,br_in(ix,j)
       enddo 
!#! Fine ciclo in js
       enddo
!$OMP END PARALLEL DO
 
       DO IX=0,LX-1 
       DSUM2(IX)=(SUM1(IX+1)-SUM1(IX)-SUM0(IX+1)+SUM0(IX))
     .          *ZS22(IX)*0.5d0
       
       DSUM3(IX) =((SUM2(IX+1)-SUM8(IX+1))*X(IX+1)-(SUM2(IX)-SUM8(IX))  
     .          *X(IX))*ZS22(IX)/X2(IX)*0.5d0
!        write(*,*) 'DF-3',ix,j,dsum2(ix),dsum3(ix)
        
       enddo

       J = JANZ(MM(IM),N)   
       J1 = JANZ(-MM(IM),-N)
       CFAK1 = N/RR*I       
!       write(*,*) 'DIAG_farad0',ix,d_t,CFAK1,DSUM2(lx-2),dsum3(lx-2)
!     .        ,vrlx,etaa
!#! Marco, 31 marzo 2017: ho inserito questo if perché c'è una
!differenza tra la prima applicazione della legge di Faraday e la
!seconda
       if (sseq .eq. 1) then
        DO IX=0,LX-1         
         CFAK = M/X2(IX)*I    
!         write(*,*) 'DF-1',ix,j,bz_in(ix,j),bt_in(ix,j)
         bt_out(IX,J) = bt_in(IX,J)+d_t*(+CFAK1*SUM3(IX)-DSUM2(IX))
         bz_out(IX,J) = bz_in(IX,J)+d_t*(+DSUM3(IX)-CFAK*SUM3(IX))
!         write(*,*) 'DF1',ix,j,bz_out(ix,j),bt_out(ix,j)
        enddo
       else if (sseq .eq. 2) then
        DO IX=0,LX-1         
         CFAK = M/X2(IX)*I    
!         write(*,*) 'DF-2',ix,j,bz_out(ix,j),bt_out(ix,j)
         bt_out(IX,J) = bt_out(IX,J)+d_t*(+CFAK1*SUM3(IX)-DSUM2(IX))
         bz_out(IX,J) = bz_out(IX,J)+d_t*(+DSUM3(IX)-CFAK*SUM3(IX))
!         write(*,*) 'DF2',ix,j,bz_out(ix,j),bt_out(ix,j)
        enddo
       endif
       ! Daniele, 18/08/2009
       ! caso generale per VRLX =/ 0 e BR =/0
       bz_out(LX,J) = bz_out(LX-1,J)*
     .        (2.*ETAA(LX)+VRLX/LX)/
     .        (2.*ETAA(LX)-VRLX/LX)
     .        + I*ETAA(LX)*N/RR*br_in(LX,J)*
     .        (2./LX)/(2.*ETAA(LX)-VRLX/LX)
       
       ! Daniele, 18/08/2009
       ! caso generale per VRLX =/ 0 e BR =/0
       bt_out(LX,J) = bt_out(LX-1,J)*  
     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
     .        + I*ETAA(LX)*M*br_in(LX,J)*
     .        (2./LX)/(2.*X2(LX)*ETAA(LX)-VRLX/LX)
       
       bt_out(-1,J) = -bt_out(0,J)
       bz_out(-1,J) = -bz_out(0,J)
       IF (J.eq.0) then
       ! Daniele, 25/11/2008
       ! caso generale per VRLX =/ 0
        bt_out(LX,J) = bt_out(LX-1,J)*
     .        (2.*X2(LX-1)*ETAA(LX) + VRLX/LX)/
     .        (2.*X2(LX)*ETAA(LX) - VRLX/LX)
     .        + E0*(2./LX)/(2.*X2(LX)*ETAA(LX)-VRLX/LX)
       ! Daniele, 25/11/2008
       ! caso generale per VRLX =/ 0
        bz_out(LX,J) = bz_out(LX-1,J)*
     .        (2.*ETAA(LX)+VRLX/LX)/
     .        (2.*ETAA(LX)-VRLX/LX)

       endif
       IF (M .eq. 0) then
         bz_out(-1,J) = bz_out(0,J)
       endif
       J1 = JANZ(-M,-N)

       call apply_solen(br_out(:,j),bt_out(:,j),bz_out(:,j),m,n)
       
       IF (J.ne.0) then
        DO IX=0,LX                                                    
         br_out(IX,J1) = CONJG(br_out(IX,J))                                    
        enddo                                                          
        DO IX=-1,LX                      
         bt_out(IX,J1) = CONJG(bt_out(IX,J))
         bz_out(IX,J1) = CONJG(bz_out(IX,J))
        enddo
       endif
       
       enddo !#! fine ciclo in N
       enddo !#! fine ciclo in M
!#! Marco, fine del ciclo m,n per prima legge di Faraday

      end subroutine faraday


      subroutine j_bestimmung

! Local variables
       integer(wp) :: im, m, n, j, j1, ix

       DO IM=0,MY        
       M = MM(IM)
       DO N=NA1(M),NE(M)
       FAK1 = N/RR          
       J = JANZ(M,N)        
       J1 = JANZ(-M,-N)     
       DO IX=0,LX-1     
        JR(IX,J) = I*(M/X2(IX)*BZN(IX,J)-FAK1*BTN(IX,J))                 
       enddo
       !#! Marco, 02/10/2014, modifica fatta con Daniele per calcolare la
       !corrente anche a LX
       !       DO 409 IX=1,LX-1                                                 
       DO IX=1,LX
        JT(IX,J) = I*FAK1*BRN(IX,J)-2.0*ZS2(IX)*(BZN(IX,J)-BZN(IX-1,J))  
        JZ(IX,J) = 2.0*ZS2(IX)/X(IX)*(X2(IX)*BTN(IX,J)-X2(IX-1)          
     .            *BTN(IX-1,J))-I*M/X(IX)*BRN(IX,J)
       enddo
       !       JT(LX,J) = (0.0,0.0)                                             
       !       JZ(LX,J) = (0.0,0.0)                                             
       JT(0,J) =(0.0,0.0)                                               
       JZ(0,J) =(0.0,0.0)                                               
       IF(M.eq.0) then
       !       JT(LX,J) = -2.*ZS2(LX)*(BZN(LX,J) - BZN(LX-1,J))                 
       !       JZ(LX,J) =  2.*ZS2(LX)*((2. - X2(LX-1))*BTN(LX,J) -              
       !     &                                          X2(LX-1)*BTN(LX-1,J))   
       JZ(0,J) = 2.*BTN(0,J)/X2(0)                                     
       endif
       IF(M.eq.1) then
        JT(0,J) = JT(1,J)
       endif
       IF(J.ne.0) then
        DO IX=0,LX                                                   
         JR(IX,J1) = CONJG(JR(IX,J))                                      
        enddo
        DO IX=0,LX                                                   
         JT(IX,J1) = CONJG(JT(IX,J))                                      
         JZ(IX,J1) = CONJG(JZ(IX,J))                                      
        enddo
       endif
       enddo
       enddo

      end subroutine j_bestimmung

      subroutine apply_solen(abrn,abtn,abzn,m,n)

       complex(wp) , intent(inout) :: abrn(0:lx)
       complex(wp) , intent(inout) :: abtn(-1:lx)
       complex(wp) , intent(inout) :: abzn(-1:lx)
       integer(wp) , intent(in) :: m, n

! Local variables

!  BR AUS BT +BZ ........................                               
       abrn(0) =(0.0,0.0)                                              
!#! Marco, calcolo anche il punto LX
       do IX=1,LX-1                                                 
       abrn(IX)=(I/ZS22(IX-1)*(-abtn(IX-1)*M-abzn(IX-1)*N*X2(IX-1)   
     .           /RR)+abrn(IX-1)*X(IX-1))/X(IX)
       enddo
       IF(M.eq.1) then
        abrn(0) = (4.0*abrn(1)-abrn(2))/3.d0
        abtn(-1) = -abtn(0)+2.d0*i*abrn(0)
       endif

      end subroutine apply_solen

      end module equations_solution


