pro ireadeq,path

;  path='./'
;#; main declarations
  zeit = double(0)
  ly = long(0)
  izz = long(0)

  abz = dcomplex(0)
  abt = dcomplex(0)
  abr = dcomplex(0)
  avz = dcomplex(0)
  avt = dcomplex(0)
  avr = dcomplex(0)

  print,'path ',path
  nsec = find_nsectionsdat(path)
  xxx = read_settings(path,nsec)

  ntp = xxx.cumulated(nsec-1)


  for k=0,nsec-1 do begin

   filein='specyl_'+strtrim(string(k+1,format='(i0)'))+'_all.dat'
   file = path+'dat/'+filein
   print,'ireadsc_file=',file

;#; opening just to read the main simulations parameters
   openr, iunit, file,/get_lun,/f77_unformatted
   print,'**************************'
   print,'opening readsc ',file
   print,'**************************'
   readu, iunit, zeit,ly,izz
   print,'izz= ',izz
;   print,'first_read',zeit,ly,izz
;#; checks
   if (ly eq 0 or izz eq 0) then begin
    print, 'ERROR! lx=0 or izz=0'
    exit
   endif
   free_lun, iunit

   lx = ly
   iz = izz

   t = make_array(ntp,/double,value=0.) 
   bz = make_array(ntp,lx+2,/dcomplex,value=0.)
   bt = make_array(ntp,lx+2,/dcomplex,value=0.)
   br = make_array(ntp,lx+1,/dcomplex,value=0.)
   vz = make_array(ntp,lx+2,/dcomplex,value=0.)
   vt = make_array(ntp,lx+2,/dcomplex,value=0.)
   vr = make_array(ntp,lx+1,/dcomplex,value=0.)

   openr, iunit, file,/get_lun,/f77_unformatted

   itp = 0
   while  ~(eof(iunit)) do begin
     readu, iunit,zeit,ly,izz
     if (itp mod 10 eq 0) then begin
      print, itp,zeit,ly,izz,ntp
     endif
     t[itp] = zeit
;     print,'second_read',izz
     for j=0,izz do begin
;      print,itp,j, avt, avz, abt, abz
      readu,iunit,avt, avz, abt, abz

      if (j eq 0) then begin
       vt(itp,0) = avt
       vz(itp,0) = avz
       bt(itp,0) = abt
       bz(itp,0) = abz
      endif

      for ix=0,ly do begin
       readu, iunit,avr, avt, avz, abr, abt, abz
       if (j eq 0) then begin
        vr(itp,ix) = avr
        vt(itp,ix+1) = avt
        vz(itp,ix+1) = avz
        br(itp,ix) = abr
        bt(itp,ix+1) = abt
        bz(itp,ix+1) = abz
       endif
      endfor
     endfor
     itp = itp+1
     if (itp ge ntp) then begin
      print,'ATTENTION! Increase the value of ntp in ireasc.pro!'
      break
     endif
   endwhile
   free_lun,iunit

   if (k eq 0) then begin
       print,'k=',k,' itp=',itp
       pebz = temporary(bz[0:itp-1,*])
       pebt = temporary(bt[0:itp-1,*])
       pebr = temporary(br[0:itp-1,*])
       pevz = temporary(vz[0:itp-1,*])
       pevt = temporary(vt[0:itp-1,*])
       pevr = temporary(vr[0:itp-1,*])
       pet = temporary(t[0:itp-1,*])
       print,'memory1',memory(/current)
   endif else begin

       print,'k=',k,' itp=',itp
       pebz=[pebz,temporary(bz[1:itp-1,*])]
       pebt=[pebt,temporary(bt[1:itp-1,*])]
       pebr=[pebr,temporary(br[1:itp-1,*])]
       pevz=[pevz,temporary(vz[1:itp-1,*])]
       pevt=[pevt,temporary(vt[1:itp-1,*])]
       pevr=[pevr,temporary(vr[1:itp-1,*])]
       pet = [pet,temporary(t[1:itp-1])]
       print,'memory2',memory(/current)

   endelse

;#; end of the cycle on the number of sections
   endfor
    eqbr = real_part(pebr)
    eqbt = real_part(pebt)
    eqbz = real_part(pebz)
    eqvr = real_part(pevr)
    eqvt = real_part(pevt)
    eqvz = real_part(pevz)
    eqt = pet

   print,'Sto salvando il file.. attendere.. ..ireadeq.pro'
   fileout = 'bv_eq.sav'
   save_file = path+'dat/'+fileout
   save,eqt,lx,eqbz,eqbt,eqbr,eqvz,eqvt,eqvr,filename=save_file

end
