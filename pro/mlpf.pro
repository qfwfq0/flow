function mlpf, mna,mnb
;#; function that computes the multiplying factor needed while doing nonlinear convolutions

   if (mna[0]*mnb[0] eq 0) then begin
     if (mna[1]*mnb[1] eq 0) then begin
       if (mna[0] eq 0 and mnb[0] eq 0) then begin
         if (mna[1] eq 0 and mnb[1] eq 0) then begin
           fact = 1.d
         endif else begin
           fact = 2.d
         endelse
       endif else begin
         fact = 2.
       endelse
     endif else begin
       if (mna[0] eq 0 and mnb[0] eq 0) then begin
         fact = 2.d
       endif else begin
         fact = 1.d
       endelse
     endelse
   endif else begin
    fact = 1.d
   endelse
    return,fact
         
end 
