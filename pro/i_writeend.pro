pro i_writeend,path,d3

;  path='./'
;  nsec=2
;  d3 = 1
  print,' d3 ',d3
  print,'path ',path

  nsectionsdat = find_nsectionsdat(path)
  print,'nsectionsdat=',nsectionsdat
  xxx = read_settings(path,nsectionsdat)
  end_time = xxx.max

  q = min(where(xxx.cumulated ge xxx.max)) + 1
  print,'apro il file numero ',q
;#; Marco, 1 marzo 2019
  if (q eq 1) then begin
   tot_itp = 0
  endif else begin
   tot_itp = xxx.cumulated(q-2)
  endelse
  print,'Parto dal tempo',tot_itp
 
  end_file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_end.dat'
  openw,end_unit,end_file,/get_lun,/f77_unformatted
;  tot_itp=0l

;  for q=1,nsectionsdat do begin

;   ntp = xxx.max(q - 1)
   file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.dat'
   print,'ireadsc_file=',file

   zeit = double(0)
   ly = 0ll
   izz = 0ll

   abz = dcomplex(0)
   abt = dcomplex(0)
   abr = dcomplex(0)
   avz = dcomplex(0)
   avt = dcomplex(0)
   avr = dcomplex(0)

;#; opening just to read the main simulations parameters
   openr, iunit, file,/get_lun,/f77_unformatted
   print,'**************************'
   print,'opening readsc ',file
   print,'**************************'
   readu, iunit, zeit,ly,izz
   print,'izz= ',zeit,ly,izz
;#; checks
   if (ly eq 0 or izz eq 0) then begin
    print, 'ERROR! lx=0 or izz=0'
;    exit
   endif
   free_lun, iunit

   lx = ly
   iz = izz
;#; Marco, mesh creation
   r1 = dindgen(lx+1)/(lx)
   r2 = dindgen(lx+2)/(lx) - 1.d/(2.d*lx)

;   t = make_array(ntp+1,/double,value=0.) 
   bz = make_array(lx+2,iz+1,/dcomplex,value=0.)
   bt = make_array(lx+2,iz+1,/dcomplex,value=0.)
   br = make_array(lx+1,iz+1,/dcomplex,value=0.)
   vz = make_array(lx+2,iz+1,/dcomplex,value=0.)
   vt = make_array(lx+2,iz+1,/dcomplex,value=0.)
   vr = make_array(lx+1,iz+1,/dcomplex,value=0.)

   openr, iunit, file,/get_lun,/f77_unformatted

   itp = 0l
   print,xxx.max
   while (itp le xxx.max ) do begin
     readu, iunit, zeit, ly,izz
     if (itp mod 10 eq 0) then begin
      print,'aa',tot_itp,zeit,ly,izz
     endif
     if (tot_itp eq xxx.max) then begin
      writeu,end_unit,zeit,long64(ly),long64(izz)
      print, 'sto scrivendo il tempo=',zeit
     endif
     for j=0,izz do begin
      readu,iunit,avt, avz, abt, abz
      if (tot_itp eq xxx.max) then begin
       writeu,end_unit,avt,avz,abt,abz
      endif

      for ix=0,ly do begin
       readu, iunit,avr, avt, avz, abr, abt, abz
       if (tot_itp eq xxx.max) then begin
        writeu,end_unit,avr,avt,avz,abr,abt,abz
       endif

      endfor
     endfor
     itp = itp + 1
     tot_itp = tot_itp + 1
     if (tot_itp ge xxx.max + 1 ) then begin
      print,'trovato!',xxx.max,' esco'
      break
     endif
   endwhile
;
;   if (tot_itp gt xxx.max + 1) then begin
;    print,'trovato! esco2'
;    break
;   endif

;  endfor
  free_lun,iunit
  free_lun,end_unit

end
