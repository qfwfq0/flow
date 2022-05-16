pro i_changeend_v2,path,d3

;  path='./'
;  nsec=2
;  d3 = 1
  print,' d3 ',d3
  print,'path ',path

  nsectionsdat = find_nsectionsdat(path)
  print,'nsectionsdat=',nsectionsdat
  xxx = read_settings(path,nsectionsdat)
  end_time = xxx.max
  print,'end_time= ',end_time
  
  end_file = path+'dat/specyl_end_changed_505.dat'
  new_izz = 505
  end_file = path+'dat/specyl_end_changed_225.dat'
  new_izz = 225
  openw,end_unit,end_file,/get_lun,/f77_unformatted
  tot_itp=0l

;#  for q=1,nsectionsdat do begin
  for q=4,4 do begin

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
    exit
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
;#;Marco, specifico
   xxx.max=700
   print,'max_t= ',xxx.max
   print,'new_izz= ',new_izz
   itp = 0l
   while (itp le xxx.max ) do begin
     readu, iunit, zeit, ly,izz
     if (itp mod 10 eq 0) then begin
      print, tot_itp,zeit,ly,izz;,ntp
     endif
     if (tot_itp eq xxx.max) then begin
      writeu,end_unit,zeit,long64(ly),long64(new_izz)
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
     for j=izz+1,new_izz do begin
      if (tot_itp eq xxx.max) then begin
       writeu,end_unit,dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d)
      endif
      for ix=0,ly do begin
       if (tot_itp eq xxx.max) then begin
        writeu,end_unit,dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d)
       endif
      endfor
     endfor


     itp = itp + 1
     tot_itp = tot_itp + 1
     if (tot_itp gt xxx.max + 1) then begin
      print,'trovato! esco'
      break
     endif
   endwhile
     if (tot_itp gt xxx.max + 1) then begin
      print,'trovato! esco2'
      break
     endif

  endfor
   free_lun,iunit
   free_lun,end_unit


end
