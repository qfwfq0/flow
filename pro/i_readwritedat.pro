pro ireadwritedat,path,nsec,mf,d3,deltat


  print,'mf ',mf,' nsec ',nsec,' d3 ',d3
  print,'path ',path

  nsectionsdat = find_nsectionsdat(path)
  print,'nsectionsdat=',nsectionsdat
  xxx = read_settings(path,nsectionsdat)

  nsec = 0

   zeit = double(0)
   ly = long(0)
   izz = long(0)

   abz = dcomplex(0)
   abt = dcomplex(0)
   abr = dcomplex(0)
   avz = dcomplex(0)
   avt = dcomplex(0)
   avr = dcomplex(0)

;#; opening just to read the main simulations parameters
   q = 1
   file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.dat'
   openr, iunit, file,/get_lun,/f77_unformatted
   print,'**************************'
   print,'opening readsc ',file
   print,'**************************'
   readu, iunit, zeit,ly,izz
   print,'izz= ',zeit,izz
;#; checks
   if (ly eq 0 or izz eq 0) then begin
    print, 'ERROR! lx=0 or izz=0'
    exit
   endif
   free_lun, iunit

   lx = ly
;#; Marco, mesh creation
   r1 = dindgen(lx+1)/(lx)
   r2 = dindgen(lx+2)/(lx) - 1.d/(2.d*lx)
;#; I change this for memory reason. I save m<=1.
   if (d3 eq 1) then begin
    iz = 90
    print,'ATTENTION! Saving only m<=1' 
   endif else begin
    iz = 10
   endelse

   deltat = 0
   t = make_array(deltat+1,/double,value=0.) 
   bz = make_array(deltat+1,lx+2,iz+1,/dcomplex,value=0.)
   bt = make_array(deltat+1,lx+2,iz+1,/dcomplex,value=0.)
   br = make_array(deltat+1,lx+1,iz+1,/dcomplex,value=0.)
   vz = make_array(deltat+1,lx+2,iz+1,/dcomplex,value=0.)
   vt = make_array(deltat+1,lx+2,iz+1,/dcomplex,value=0.)
   vr = make_array(deltat+1,lx+1,iz+1,/dcomplex,value=0.)

   openr, iunit, file,/get_lun,/f77_unformatted

   itp = 0
   ndelta = 0
  readu, iunit,zeit,ly,izz
     if (pitp mod 10 eq 0) then begin
      print, itp,zeit,ly,izz,ntp
     endif
     t[pitp] = zeit
;     print,'second_read',izz
     for j=0,izz do begin
;      print,itp,j, avt, avz, abt, abz
      readu,iunit,avt, avz, abt, abz

      if (j lt iz) then begin
      vt(pitp,0,j) = avt
      vz(pitp,0,j) = avz
      bt(pitp,0,j) = abt
      bz(pitp,0,j) = abz
      endif

;      if (j eq 0) then begin
;       vt0(pitp,0) = avt
;       vz0(itp,0) = avz
;       bt0(itp,0) = abt
;       bz0(itp,0) = abz
;      endif

      for ix=0,ly do begin
       readu, iunit,avr, avt, avz, abr, abt, abz

       if (j lt iz) then begin
       vr(pitp,ix,j) = avr
       vt(pitp,ix+1,j) = avt
       vz(pitp,ix+1,j) = avz
       br(pitp,ix,j) = abr
       bt(pitp,ix+1,j) = abt
       bz(pitp,ix+1,j) = abz
       endif

;       if (j eq 0) then begin
;        vr0(itp,ix) = avr
;        vt0(itp,ix+1) = avt
;        vz0(itp,ix+1) = avz
;        br0(itp,ix) = abr
;        bt0(itp,ix+1) = abt
;        bz0(itp,ix+1) = abz
;       endif

      endfor

     free_lun,iunit
     print,pitp

     print,'Sto salvando il file.. attendere.. ..ireadsc.pro'
     filew = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_equil.dat'
     openw, iunitw, filew,/get_lun,/f77_unformatted

     writeu, iunit,zeit
     writeu, iunit,ly
     writeu, iunit,izz
     for k = 0, 101 do begin
     endfor

     free_lun,iunitw

end
