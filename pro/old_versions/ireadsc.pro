pro ireadsc,path,fileout,filein,sec,nsec,mf,d3

;   path='./'
;#; main declarations
   print,'mfi ',mf,' nsec',nsec,' d3',d3
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
;   set_file=path+'for/settings.blc.for'
;   openr,lun,set_file,/get_lun
;   char=make_array(8,/string)
;   readf,lun,char
;;#; extract, from settings.blc.for, the maximum time that needs to be analyzed
;   dum=char(5)
;   a=strpos(dum,'/')
;   dumm=strmid(dum,a+1)
;   a=strpos(dumm,'/')
;   dummm=strmid(dumm,0,a)
;;#; ntp is defined from input data
;   ntp=long(dummm)+1
   nsectionsdat = find_nsectionsdat(path)
   print,'nsectionsdat=',nsectionsdat
   xxx = read_settings(path,nsectionsdat)

   if (nsec gt nsectionsdat) then begin
    print,'tryying to compute a file that does not exist yet! ',nsec,nsectionsdat
    stop
   endif
   ntp = xxx.computed(nsec-1)
   file = path+'dat/'+filein
   print,'ireadsc_file=',file

;#; opening just to read the main simulations parameters
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
;#; I change this for memory reason. I save m<=1.
       print,' ' 
       print,'d3d3d3',d3
   print,'ATTENTION! Saving only m<=1' 
   if (d3 eq 1) then begin
    iz = 90
   endif else begin
    iz = 10
   endelse

   t = make_array(ntp+1,/double,value=0.) 
   bz = make_array(ntp+1,lx+2,iz+1,/dcomplex,value=0.)
   bt = make_array(ntp+1,lx+2,iz+1,/dcomplex,value=0.)
   br = make_array(ntp+1,lx+1,iz+1,/dcomplex,value=0.)
   vz = make_array(ntp+1,lx+2,iz+1,/dcomplex,value=0.)
   vt = make_array(ntp+1,lx+2,iz+1,/dcomplex,value=0.)
   vr = make_array(ntp+1,lx+1,iz+1,/dcomplex,value=0.)

;#; arrays to save the equilibrium part of the field.
   bz0 = make_array(ntp+1,lx+2,/dcomplex,value=0.)
   bt0 = make_array(ntp+1,lx+2,/dcomplex,value=0.)
   br0 = make_array(ntp+1,lx+1,/dcomplex,value=0.)
   vz0 = make_array(ntp+1,lx+2,/dcomplex,value=0.)
   vt0 = make_array(ntp+1,lx+2,/dcomplex,value=0.)
   vr0 = make_array(ntp+1,lx+1,/dcomplex,value=0.)

   openr, iunit, file,/get_lun,/f77_unformatted

   itp = 0
   while ~(eof(iunit)) do begin
     readu, iunit,zeit,ly,izz
;     if (itp mod 10 eq 0) then begin
      print, itp,zeit,ly,izz,ntp
;     endif
     t[itp] = zeit
;     print,'second_read',izz
     for j=0,izz do begin
;      print,itp,j, avt, avz, abt, abz
      readu,iunit,avt, avz, abt, abz

      if (j lt iz) then begin
      vt(itp,0,j) = avt
      vz(itp,0,j) = avz
      bt(itp,0,j) = abt
      bz(itp,0,j) = abz
      endif

      if (j eq 0) then begin
       vt0(itp,0) = avt
       vz0(itp,0) = avz
       bt0(itp,0) = abt
       bz0(itp,0) = abz
      endif

      for ix=0,ly do begin
       readu, iunit,avr, avt, avz, abr, abt, abz

       if (j lt iz) then begin
       vr(itp,ix,j) = avr
       vt(itp,ix+1,j) = avt
       vz(itp,ix+1,j) = avz
       br(itp,ix,j) = abr
       bt(itp,ix+1,j) = abt
       bz(itp,ix+1,j) = abz
       endif

       if (j eq 0) then begin
        vr0(itp,ix) = avr
        vt0(itp,ix+1) = avt
        vz0(itp,ix+1) = avz
        br0(itp,ix) = abr
        bt0(itp,ix+1) = abt
        bz0(itp,ix+1) = abz
       endif

      endfor
     endfor
     if (itp ge ntp) then begin
      print,'ATTENTION! Increase the value of ntp in ireasc.pro!'
      break
     endif
     itp = itp+1
   endwhile
   free_lun,iunit

   t = temporary(t[0:itp-1])
   bz = temporary(bz[0:itp-1,*,*])
   bt = temporary(bt[0:itp-1,*,*])
   br = temporary(br[0:itp-1,*,*])
   vz = temporary(vz[0:itp-1,*,*])
   vt = temporary(vt[0:itp-1,*,*])
   vr = temporary(vr[0:itp-1,*,*])

   bz0 = temporary(bz0[0:itp-1,*])
   bt0 = temporary(bt0[0:itp-1,*])
   br0 = temporary(br0[0:itp-1,*])
   vz0 = temporary(vz0[0:itp-1,*])
   vt0 = temporary(vt0[0:itp-1,*])
   vr0 = temporary(vr0[0:itp-1,*])

   print,'Sto salvando il file.. attendere.. ..ireadsc.pro'
   save_file = path+'dat//'+fileout
   save,t,lx,bz,bt,br,vz,vt,vr,filename=save_file

;   print,'Sto salvando il file.. attendere.. ..ireadsc.pro'
;   fileout = 'bv_eq'+strtrim(string(nsec,format='(i0)'))+'.sav'
;   save_file = path+'dat/bv_eq_'+string(nsec,format='(i0)')+'.sav'
;   save,t,lx,bz0,bt0,br0,vz0,vt0,vr0,filename=save_file

;#; save the field in modulus&phase notation
   if (mf eq 1) then begin
    bfield_mf_out=path+'dat/imf_bprofiles_'+string(sec,format='(i0)')+'.sav'
    call_procedure,'mf1',t,br,bt,bz,bfield_mf_out

    vfield_mf_out=path+'dat/imf_vprofiles_'+string(sec,format='(i0)')+'.sav'
    call_procedure,'mf1',t,vr,vt,vz,vfield_mf_out
   endif

end
