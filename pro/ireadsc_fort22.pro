pro ireadsc_fort22,path,nsec,mf,d3

  print,'mf ',mf,' nsec ',nsec,' d3 ',d3
  print,'path ',path

  nsectionsdat = find_nsectionsdat(path)
  print,'nsectionsdat=',nsectionsdat
  xxx = read_settings(path,nsectionsdat)

  if (nsec gt nsectionsdat) then begin
   print,'tryying to compute a file that does not exist yet! ',nsec,nsectionsdat
   stop
  endif


  for q=2,2 do begin

   compute_it=test_up_to_date_files(path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav',path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_idl.dat')
   print,'ireadsc_fort',q,compute_it
;   if (compute_it gt 0) then begin
   ntp = xxx.computed(q - 1)
   file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_idl.dat'
   file = path+'dat/specyl_21_idl.dat'
   print,'ireadsc_file=',file

   zeit = double(0)
   ly = 0ll
   izz = 0ll
   max_time = 0ll
   time = 0ll

;#; opening just to read the main simulations parameters
   openr, iunit, file,/get_lun,/f77_unformatted
   print,'**************************'
   print,'opening readsc ',file
   print,'**************************'
   readu, iunit, ly
   readu, iunit, izz
   readu, iunit, max_time
   readu, iunit, time
   print,'izz= ',ly,izz,max_time,time
;#; checks
   if (ly eq 0 or izz eq 0) then begin
    print, 'ERROR! lx=0 or izz=0'
    exit
   endif

   time=5000
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
    iz = izz
   endelse

;   t = make_array(max_time+1,/double,value=0.) 
;   bz = make_array(max_time+1,lx+2,iz+1,/dcomplex,value=0.)
;   bt = make_array(max_time+1,lx+2,iz+1,/dcomplex,value=0.)
;   br = make_array(max_time+1,lx+1,iz+1,/dcomplex,value=0.)
;   vz = make_array(max_time+1,lx+2,iz+1,/dcomplex,value=0.)
;   vt = make_array(max_time+1,lx+2,iz+1,/dcomplex,value=0.)
;   vr = make_array(max_time+1,lx+1,iz+1,/dcomplex,value=0.)
   t = make_array(time+1,/double,value=0.) 
   bz = make_array(time+1,lx+2,iz+1,/dcomplex,value=0.)
   bt = make_array(time+1,lx+2,iz+1,/dcomplex,value=0.)
   br = make_array(time+1,lx+1,iz+1,/dcomplex,value=0.)
   vz = make_array(time+1,lx+2,iz+1,/dcomplex,value=0.)
   vt = make_array(time+1,lx+2,iz+1,/dcomplex,value=0.)
   vr = make_array(time+1,lx+1,iz+1,/dcomplex,value=0.)
   readu, iunit, t
   readu, iunit, br
   readu, iunit, bt
   readu, iunit, bz
   readu, iunit, vr
   readu, iunit, vt
   readu, iunit, vz

   t = t(0:time-1)
   br = br(0:time-1,*,*)
   bt = bt(0:time-1,*,*)
   bz = bz(0:time-1,*,*)
   vr = vr(0:time-1,*,*)
   vt = vt(0:time-1,*,*)
   vz = vz(0:time-1,*,*)

   print,'Sto salvando il file.. attendere.. ..ireadsc.pro'
   save_file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav'
   save,t,lx,r1,r2,bz,bt,br,vz,vt,vr,filename=save_file

   spawn, "rm -f "+file
;   print,bt[0,101,1]
;   print,'****'
;   print,'Sto salvando il file.. attendere.. ..ireadsc.pro'
;   fileout = 'bv_eq'+strtrim(string(nsec,format='(i0)'))+'.sav'
;   save_file = path+'dat/bv_eq_'+string(nsec,format='(i0)')+'.sav'
;   save,t,lx,bz0,bt0,br0,vz0,vt0,vr0,filename=save_file

;#; remove the specyl_?_idl.dat file
;    spawn, 'rm -f '+strtrim(file)


stop
;#; save the field in modulus&phase notation
   print,'test',mf
   if (mf eq 1) then begin
    bfield_mf_out=path+'dat/imf_bprofiles_'+string(q,format='(i0)')+'.sav'
    call_procedure,'mf1',t,br,bt,bz,bfield_mf_out

    vfield_mf_out=path+'dat/imf_vprofiles_'+string(q,format='(i0)')+'.sav'
    call_procedure,'mf1',t,vr,vt,vz,vfield_mf_out
   endif

;#; end of the compute_it if
;   endif

;#; Marco, treat the case in which the imf_?profiles_?.sav do not exist
   if (mf eq 1) then begin
    bfield_mf_out=path+'dat/imf_bprofiles_'+string(q,format='(i0)')+'.sav'
    vfield_mf_out=path+'dat/imf_vprofiles_'+string(q,format='(i0)')+'.sav'
    info_ex = file_info(bfield_mf_out)
    print,'CHECK_existence_ireadsc_fort',info_ex.exists
    if (info_ex.exists eq 0)then begin
     compute_merge=test_up_to_date_files(path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.dat',path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav')
     if (compute_merge eq 1) then begin
      ;#; nothing to do
     endif else begin
      save_file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav'
      restore,save_file
      call_procedure,'mf1',t,br,bt,bz,bfield_mf_out
      call_procedure,'mf1',t,vr,vt,vz,vfield_mf_out
     endelse
    endif
;#; Marco, controllo il file di velocità
     info_exv = file_info(vfield_mf_out)
     
     if (info_exv.exists eq 0)then begin
      save_file = path+'dat/specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav'
      restore,save_file
      call_procedure,'mf1',t,vr,vt,vz,vfield_mf_out
     endif
   endif

  endfor

;#; 
  compute_merge=test_up_to_date_files(path+'dat/specyl_'+strtrim(string(nsec,format='(i0)'))+'_all.sav',path+'dat/imf_bprofiles.sav')
  if (compute_merge eq 1) then begin
   print,'already merged b'
  endif else begin
   call_procedure,'imergebmf',path,d3,mf
  endelse

  compute_merge=test_up_to_date_files(path+'dat/specyl_'+strtrim(string(nsec,format='(i0)'))+'_all.sav',path+'dat/imf_vprofiles.sav')
  if (compute_merge eq 1) then begin
   print,'already merged v'
  endif else begin
   call_procedure,'imergevmf',path,d3,mf
  endelse


end
