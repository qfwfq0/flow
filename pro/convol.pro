
function convol mcap,ncap,my,mm,nanf,nz

path='./'
;input='spectrum.sav'
;restore,path+'dat/'+strtrim(input)

;mcap = 3
;ncap = -19

;#; consistency check.
if (mcap gt my) then begin
 print,'wrong m',mcap
 stop
endif
if (ncap lt nanf(mcap) or ncap gt nanf(mcap)+nz(mcap)-1) then begin
 print,'wrong n',ncap
 stop
endif
emme = make_array(2,my+1,/integer,value=-1000)
enne = make_array(2,max(nz),/integer,value=-1000)
c1=0
c2=0
c3=0
c4=0
;#; parte mcap=0
if (mcap eq 0) then begin
   mtilda = 0
   emme(0,mtilda) = 0
   emme(1,mtilda) = 0
 ;#; part where i find ntilda+enne = ncap
   for ntilda = nanf(0),nanf(0) + nz(0)-1 do begin
     ntemp = ncap - ntilda
     if (ntemp lt nanf(0) or ntemp gt nanf(0)+nz(0)-1) then begin
      continue
     endif else begin
      print,'n 0',ntilda,ntemp
      enne(0,c1) = ntilda
      enne(1,c1) = ntemp
      c1 = c1 + 1
     endelse
   endfor

;#; parte mcap,ncap /= 0
endif else begin

 for mtilda=0,mm(my) do begin
;#; case when mtilda eq 0
   if (mtilda eq 0) then begin
    emme(0,mtilda) = 0
    emme(1,mtilda) = mcap - 0 
;#; part where i find ntilda+enne = ncap
     for ntilda = nanf(0),nanf(0) + nz(0)-1 do begin
      if (ntilda ge 0) then begin ;#; special case because I know that I'm dealing with m=0
       continue
      endif else begin
       ntemp = ncap - ntilda
       if (ntemp lt nanf(emme(mtilda)) or ntemp gt nanf(emme(mtilda))+nz(emme(mtilda))-1) then begin
        continue
       endif else begin
        print,'mtilda, m= 0,',emme(mtilda),' n,ntilda=',ntilda,ntemp
        enne(0,c2) = ntilda
        enne(1,c2) = ntemp
        c2 = c2 + 1
       endelse
      endelse
     endfor
   endif else begin
   mtemp = mcap - mtilda
   if (mtemp lt 0) then begin
    continue
   endif else begin
;#; caso mtemp greater than zero

;#; caso mtemp equal zero
    if (mtemp eq 0) then begin
     print,'cacca mtilda, mtemp',mtilda,mtemp
     emme(0,mtilda) = mtilda
     emme(1,mtilda) = mtemp
;#; part where i find ntilda+enne = ncap
     for ntilda = nanf(mtilda),nanf(mtilda) + nz(mtilda) -1 do begin
      if (ntilda eq 0) then begin 
       continue
      endif else begin
       ntemp = ncap - ntilda
       if (ntemp lt nanf(mtemp) or ntemp gt nanf(mtemp)+nz(mtemp)-1 or ntemp ge 0) then begin
        continue
       endif else begin
        enne(0,c3) = ntilda
        enne(1,c3) = ntemp
        c3 = c3 + 1
        print, 'cacca, ntilda,ntemp ',ntilda,ntemp
       endelse
      endelse
     endfor
    endif else begin
     emme(0,mtilda) = mtilda
     emme(1,mtilda) = mtemp
     print,' easy mtilda,mtemp',mtilda,mtemp
  ;#; part where i find ntilda+enne = ncap
     for ntilda = nanf(mtilda),nanf(mtilda) + nz(mtilda)-1 do begin
      ntemp = ncap - ntilda
      if (ntemp lt nanf(mtemp) or ntemp gt nanf(mtemp)+nz(mtemp)-1) then begin
       continue
      endif else begin
       print,'n easy',ntilda,ntemp
        enne(0,c4) = ntilda
        enne(1,c4) = ntemp
        c4 = c4 + 1
      endelse
     endfor
    endelse
   endelse
  endelse
 endfor
endelse


end
