pro ireadsc2,path,nsec,mf,d3

  print,'mf ',mf,' nsec ',nsec,' d3 ',d3
  print,'path ',path

  nsectionsdat = find_nsectionsdat(path)
  print,'nsectionsdat=',nsectionsdat
  xxx = read_settings(path,nsectionsdat)

  if (nsec gt nsectionsdat) then begin
   print,'tryying to compute a file that does not exist yet! ',nsec,nsectionsdat
   stop
  endif


  print,xxx
  for q=1,nsec do begin

   compute_it=up_to_date_all(path,q)
   print,'ireadsc2',q,compute_it
   if (compute_it gt 0) then begin
   ntp = xxx.computed(q - 1)
   file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.dat'
   print,'ireadsc_file=',file

   zeit = double(0)
   ly = 0ll
   izz = 0ll
;#; Marco, a volte bisogna cambiare tipo di integer
;   ly = 0l
;   izz = 0l

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
;#; Marco, 02/2017 dopo ifort17  -- falso allarme??
;   readu, iunit, zeit
;   readu, iunit, ly,izz
   print,'izz= ',zeit,ly,izz
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
    iz = izz
;    print,'ATTENTION! Saving only m<=1' 
;    iz = 135
;    print,'ATTENTION! Saving only m<=2' 
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
;#; Marco, 02/2017 dopo ifort17
;   readu, iunit, zeit,ly,izz
     readu, iunit, zeit, ly,izz
     if (itp mod 10 eq 0) then begin
      print, itp,zeit,ly,izz,ntp
     endif
     t[itp] = zeit
;     print,'second_read',izz
     for j=0,izz do begin
;      print,itp,j, avt, avz, abt, abz
      readu,iunit,avt, avz, abt, abz

      if (j le iz) then begin
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

       if (j le iz) then begin
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
   save_file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav'
   save,t,lx,r1,r2,bz,bt,br,vz,vt,vr,filename=save_file

;   print,bt[0,101,1]
;   print,'****'
;   print,'Sto salvando il file.. attendere.. ..ireadsc.pro'
;   fileout = 'bv_eq'+strtrim(string(nsec,format='(i0)'))+'.sav'
;   save_file = path+'dat/bv_eq_'+string(nsec,format='(i0)')+'.sav'
;   save,t,lx,bz0,bt0,br0,vz0,vt0,vr0,filename=save_file

;#; save the field in modulus&phase notation
   if (mf eq 1) then begin
    bfield_mf_out=path+'dat/imf_bprofiles_'+string(q,format='(i0)')+'.sav'
    call_procedure,'mf1',t,br,bt,bz,bfield_mf_out

    vfield_mf_out=path+'dat/imf_vprofiles_'+string(q,format='(i0)')+'.sav'
    call_procedure,'mf1',t,vr,vt,vz,vfield_mf_out
   endif

;#; end of the compute_it if
   endif
;#; Marco, treat the case in which the imf_?profiles_?.sav do not exist
   if (mf eq 1) then begin
    bfield_mf_out=path+'dat/imf_bprofiles_'+string(q,format='(i0)')+'.sav'
    vfield_mf_out=path+'dat/imf_vprofiles_'+string(q,format='(i0)')+'.sav'
    info_ex = file_info(bfield_mf_out)
    if (info_ex.exists eq 0)then begin
     save_file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav'
     restore,save_file
     call_procedure,'mf1',t,br,bt,bz,bfield_mf_out
     call_procedure,'mf1',t,vr,vt,vz,vfield_mf_out
    endif
   endif

  endfor
end
