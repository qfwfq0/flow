function convo,path,mcap,ncap,d3

;path='./'
;input='spectrum.sav'
;restore,path+'dat/'+strtrim(input)
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse

;mcap=0
;ncap=-7
;#; consistency check.
;print,'mcap , ncap',mcap,ncap
if (mcap gt my) then begin
 print,'wrong m',mcap
 stop
endif
if (ncap lt nanf(mcap) or ncap gt nanf(mcap)+nz(mcap)-1) then begin
 print,nanf(mcap), nz(mcap),nanf(mcap)+nz(mcap)-1
 print,'wrong n',ncap
 stop
endif

;#; 3D case

if (d3 eq 1) then begin
janzt = make_array(3,300,/integer,value=-1000)
c5=0
;#; parte mcap=0
if (mcap eq 0) then begin
   mtilda = 0
   mtemp = 0
 ;#; part where i find ntilda+enne = ncap
   for ntilda = nanf(0),nanf(0) + nz(0)-1 do begin
     ntemp = ncap - ntilda
     if (ntemp lt nanf(0) or ntemp gt nanf(0)+nz(0)-1) then begin
      continue
     endif else begin
      tildajanz = janz(mtilda,ntilda,my,mm,nanf,nz)
      tempjanz = janz(mtemp,ntemp,my,mm,nanf,nz)
      janzt(0,c5) = tildajanz
      janzt(1,c5) = tempjanz
      c5 = c5 + 1
     endelse
   endfor

;#; parte mcap,ncap /= 0
endif else begin

 for mtilda=0,mm(my) do begin
;#; case when mtilda eq 0
   if (mtilda eq 0) then begin
;#; part where i find ntilda+enne = ncap
     for ntilda = nanf(0),nanf(0) + nz(0)-1 do begin
;#; Marco, 09/12/2014: questa formulazione esclude il modo (0,0) dalla convoluzione! No!
;      if (ntilda ge 0) then begin ;#; special case because I know that I'm dealing with m=0
      if (ntilda gt 0) then begin ;#; special case because I know that I'm dealing with m=0
       continue
      endif else begin
       ntemp = ncap - ntilda
       if (ntemp lt nanf(mcap) or ntemp gt nanf(mcap)+nz(mcap)-1) then begin
        continue
       endif else begin
        tildajanz = janz(0,ntilda,my,mm,nanf,nz)
        tempjanz = janz(mcap,ntemp,my,mm,nanf,nz)
        janzt(0,c5) = tildajanz
        janzt(1,c5) = tempjanz
        c5 = c5 + 1
       endelse
      endelse
     endfor
;#; fine della parte mtida=0, in cui cerco convoluzioni con i modi m=0
   endif else begin

   mtemp = mcap - mtilda
   if (mtemp lt 0) then begin
    continue
   endif else begin
;#; caso mtemp ge than zero
;#; caso mtemp equal zero
    if (mtemp eq 0) then begin
;#; Marco, perche aveviscritto questo????
;    print,mtilda,ntilda
;      tildajanz = janz(mtilda,ntilda,my,mm,nanf,nz)
;;#; the only convolution mtemp=0 with (mcap,ncap) is ntemp=0
;      tempjanz = janz(mtemp,0,my,mm,nanf,nz)
;      janzt(0,c5) = tildajanz
;      janzt(1,c5) = tempjanz
;      c5 = c5 + 1

;#; part where i find ntilda+enne = ncap
     for ntilda = nanf(mtilda),nanf(mtilda) + nz(mtilda) -1 do begin
      if (ntilda eq 0) then begin 
;#; (mtilda,ntilda) + (mtemp,ntemp) = (mcap,ncap)
       ntemp = ncap - ntilda
       tildajanz = janz(mtilda,ntilda,my,mm,nanf,nz)
       tempjanz = janz(mtemp,ntemp,my,mm,nanf,nz)
       janzt(0,c5) = tildajanz
       janzt(1,c5) = tempjanz
       c5 = c5 + 1
      endif else begin
       ntemp = ncap - ntilda
       if (ntemp lt nanf(mtemp) or ntemp gt nanf(mtemp)+nz(mtemp)-1 or ntemp ge 0) then begin
        continue
       endif else begin
       tildajanz = janz(mtilda,ntilda,my,mm,nanf,nz)
       tempjanz = janz(mtemp,ntemp,my,mm,nanf,nz)
       janzt(0,c5) = tildajanz
       janzt(1,c5) = tempjanz
       c5 = c5 + 1
       endelse
      endelse
     endfor
;#; Marco, fine parte dove mtemp uguale a zero
    endif else begin
  ;#; part where i find ntilda+enne = ncap
     for ntilda = nanf(mtilda),nanf(mtilda) + nz(mtilda)-1 do begin
      ntemp = ncap - ntilda
      if (ntemp lt nanf(mtemp) or ntemp gt nanf(mtemp)+nz(mtemp)-1) then begin
       continue
      endif else begin
      tildajanz = janz(mtilda,ntilda,my,mm,nanf,nz)
      tempjanz = janz(mtemp,ntemp,my,mm,nanf,nz)
      janzt(0,c5) = tildajanz
      janzt(1,c5) = tempjanz
      c5 = c5 + 1
      endelse
     endfor
    endelse
   endelse
  endelse
 endfor

endelse

janzt = janzt[*,0:n_elements(where(janzt[0,*] ge -40))-1]
;sorting for dyad representing the same couple
for i = 0, n_elements(janzt[0,*])-1 do begin
 delta = 0.5 * ( max(janzt[0:1,i]) * (max(janzt[0:1,i]) + 1) ) + min(janzt[0:1,i])
 janzt[2,i] = delta
endfor
;#; find unique dyad in janzt
aaa=sort(janzt[2,*])
bbb= uniq(janzt[2,aaa])
janzf= janzt[*, aaa(bbb)]


return,janzf
endif else begin
;#; 2D case
 janzf = [0,janz(mcap,ncap,my,mm,nanf,nz)]
 return,janzf
endelse

end
