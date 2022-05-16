      module modpoinc
       
      use defs
 
      contains

         function mnum(j,nz,nanf,my,d3)
          
           implicit none
           
           integer, intent(in) :: d3, j,my
           integer, intent(in) :: nz(0:my)
           integer, intent(in) :: nanf(0:my)
           integer :: mnum(0:1), mnum2(0:1)
           integer :: sizea, mm,q,nn
           integer :: cnz(0:my)
            
           mnum2(0) = 0
           mnum2(1) = 1
           if (d3 .eq. 1) then 
!          ;#; compute the cumulated nz array, cnz
              if (j .eq. 0) then 
               mnum2(0) = 0 ; mnum2(1) = 0
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
                mnum2(0) = q
                exit
               endif
              enddo
              if (mnum2(0) .eq. 0) then
               mnum2(1) = nanf(mnum2(0)) + j - 1
              else
               mnum2(1) = nanf(mnum2(0)) +
     .          j - 1 - cnz(mnum2(0)-1) - 1
              endif
              endif
!              write(*,*) j,mnum2
            else
!          ;#; 2D case
              if (j .eq. 0) then
              mnum2(0) = 0 ; mnum2(1) = 0
              else
               mnum2(0) = j ; mnum2(1) = nanf(j)
              endif
            endif
            mnum = mnum2
          end function 

          function janz2(am,an,anz,ananf)
            
            integer(wp), intent(in) :: am,an
            integer(wp), intent(in) :: anz(:),ananf(:)
            integer(wp) :: aj,mi
!            integer(wp), intent(out) :: janz2

              aj = 0
              if (am .eq. 0 .and. an .eq. 0) then 
               janz2 = 0
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

              janz2 = aj
          end

      end module modpoinc
