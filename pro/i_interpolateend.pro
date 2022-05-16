pro i_interpolateend,path,d3

;  path='./'
;  nsec=2
;  d3 = 1
  print,' d3 ',d3
  print,'path ',path

  nsectionsdat = find_nsectionsdat(path)
  print,'nsectionsdat=',nsectionsdat
  xxx = read_settings(path,nsectionsdat)
  end_time = xxx.max


  end_file = path+'dat/specyl_end_changed_interpolated.dat'
  new_izz = 225
  new_ly = 200
  openw,end_unit,end_file,/get_lun,/f77_unformatted
  tot_itp=0l

  for q=1,nsectionsdat do begin

   ntp = xxx.max(q - 1)
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
   bz = make_array(ly+2,new_izz+1,/dcomplex,value=0.)
   bt = make_array(ly+2,new_izz+1,/dcomplex,value=0.)
   br = make_array(ly+1,new_izz+1,/dcomplex,value=0.)
   vz = make_array(ly+2,new_izz+1,/dcomplex,value=0.)
   vt = make_array(ly+2,new_izz+1,/dcomplex,value=0.)
   vr = make_array(ly+1,new_izz+1,/dcomplex,value=0.)
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
   new_bz = make_array(new_ly+2,new_izz+1,/dcomplex,value=0.)
   new_bt = make_array(new_ly+2,new_izz+1,/dcomplex,value=0.)
   new_br = make_array(new_ly+1,new_izz+1,/dcomplex,value=0.)
   new_vz = make_array(new_ly+2,new_izz+1,/dcomplex,value=0.)
   new_vt = make_array(new_ly+2,new_izz+1,/dcomplex,value=0.)
   new_vr = make_array(new_ly+1,new_izz+1,/dcomplex,value=0.)

   openr, iunit, file,/get_lun,/f77_unformatted

   itp = 0l
   while (itp le xxx.max ) do begin
     readu, iunit, zeit, ly,izz
     if (itp mod 10 eq 0) then begin
      print, tot_itp,itp,zeit,ly,izz,ntp
     endif
     if (tot_itp eq xxx.max) then begin
      writeu,end_unit,zeit,long64(ly),long64(izz)
     endif
     for j=0,izz do begin
      readu,iunit,avt, avz, abt, abz
      
      if (tot_itp eq xxx.max) then begin
        vt[0,j] = avt
        vz[0,j] = avz
        bt[0,j] = abt
        bz[0,j] = abz
;       writeu,end_unit,avt,avz,abt,abz
      endif
      for ix=0,ly do begin
       readu, iunit,avr, avt, avz, abr, abt, abz
       if (tot_itp eq xxx.max) then begin
        vr[ix,j] = avr
        br[ix,j] = abr
        vt[ix+1,j] = avt
        vz[ix+1,j] = avz
        bt[ix+1,j] = abt
        bz[ix+1,j] = abz
;        writeu,end_unit,avr,avt,avz,abr,abt,abz
       endif
      endfor
     endfor
     for j=izz+1,new_izz do begin
      if (tot_itp eq xxx.max) then begin
        vt[0,j] = dcomplex(0.d)
        vz[0,j] = dcomplex(0.d)
        bt[0,j] = dcomplex(0.d)
        bz[0,j] = dcomplex(0.d)
;       writeu,end_unit,dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d)
      endif
      for ix=0,ly do begin
       if (tot_itp eq xxx.max) then begin
        vr[ix,j] = dcomplex(0.d)
        br[ix,j] = dcomplex(0.d)
        vt[ix+1,j] = dcomplex(0.d)
        vz[ix+1,j] = dcomplex(0.d)
        bt[ix+1,j] = dcomplex(0.d)
        bz[ix+1,j] = dcomplex(0.d)
;        writeu,end_unit,dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d),dcomplex(0.d)
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

;#; interpolation part
;#; Marco, new mesh creation
   new_r1 = dindgen(new_ly+1)/(new_ly)
   new_r2 = dindgen(new_ly+2)/(new_ly) - 1.d/(2.d*new_ly)

   for j=0,new_izz do begin 
    new_br[*,j] = interpol(reform(br[*,j]),r1,new_r1) 
    new_vr[*,j] = interpol(reform(vr[*,j]),r1,new_r1) 
    new_bt[*,j] = interpol(reform(bt[*,j]),r2,new_r2) 
    new_vt[*,j] = interpol(reform(vt[*,j]),r2,new_r2) 
    new_bz[*,j] = interpol(reform(bz[*,j]),r2,new_r2) 
    new_vz[*,j] = interpol(reform(vz[*,j]),r2,new_r2) 
   endfor

   writeu,end_unit,zeit,long64(new_ly),long64(new_izz)
   for j=0,new_izz do begin
    writeu,end_unit,new_vt[0,j],new_vz[0,j],new_bt[0,j],new_bz[0,j]
    for ix=0,new_ly do begin
     writeu,end_unit,new_vr[ix,j],new_vt[ix,j],new_vz[ix,j],new_br[ix,j],new_bt[ix,j],new_bz[ix,j]
    endfor
   endfor

   free_lun,iunit
   free_lun,end_unit


end
