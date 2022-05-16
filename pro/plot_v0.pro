pro plot_v0, path, itp
;#; program to contour vt00,vz00 in two section at z=0, th=0

 xxx = read_settings(path,nsectionsdat)

 print,itp
 help,xxx,/str
 psym=0

 idl_file=path+'dat/ijprofiles.sav'
 print,'restoring: '+idl_file
 exist_file = file_test(idl_file)
 if (exist_file ne 1) then begin
  print,'need to create the ijprofiles.sav file!'
  call_procedure,'isavejprofiles',path
 endif
 restore,idl_file
 
 if (itp lt xxx.min) then begin
  print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
  stop
 endif

 ditp=fix(itp-xxx.min) 
 jzt=make_array(2,n_elements(rpjz[0,*,0]),n_elements(rpjz[0,0,*]),/double,value=0.)
 jtt=make_array(2,n_elements(rpjt[0,*,0]),n_elements(rpjt[0,0,*]),/double,value=0.)
 jrt=make_array(2,n_elements(rpjr[0,*,0]),n_elements(rpjr[0,0,*]),/double,value=0.)
 
 jzt[0,*,*]=reform(rpjz[ditp,*,*]) 
 jzt[1,*,*]=reform(ipjz[ditp,*,*]) 
 jtt[0,*,*]=reform(rpjt[ditp,*,*]) 
 jtt[1,*,*]=reform(ipjt[ditp,*,*]) 
 jrt[0,*,*]=reform(rpjr[ditp,*,*]) 
 jrt[1,*,*]=reform(ipjr[ditp,*,*]) 
 
;#; write the current field using the modulus-phase convention
 mf_jzt=make_array(2,n_elements(rpjz[0,*,0]),n_elements(rpjz[0,0,*]),/double,value=0.)
 mf_jtt=make_array(2,n_elements(rpjt[0,*,0]),n_elements(rpjt[0,0,*]),/double,value=0.)
 mf_jrt=make_array(2,n_elements(rpjr[0,*,0]),n_elements(rpjr[0,0,*]),/double,value=0.)
 mf_jzt[0,*,*] = sqrt( jzt[0,*,*]^2.d + jzt[1,*,*]^2.d )
 mf_jtt[0,*,*] = sqrt( jtt[0,*,*]^2.d + jtt[1,*,*]^2.d )
 mf_jrt[0,*,*] = sqrt( jrt[0,*,*]^2.d + jrt[1,*,*]^2.d )
 mf_jzt[1,*,*] = atan( jzt[1,*,*],jzt[0,*,*] )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_jzt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_jzt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_jzt[1,aaa[0,k],aaa[1,k]] = mf_jzt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

 mf_jtt[1,*,*] = atan( jtt[1,*,*],jtt[0,*,*] )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_jtt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_jtt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_jtt[1,aaa[0,k],aaa[1,k]] = mf_jtt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

 mf_jrt[1,*,*] = atan( jrt[1,*,*],jrt[0,*,*] )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_jrt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_jrt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_jrt[1,aaa[0,k],aaa[1,k]] = mf_jrt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

;#; saving 
 sav_file = path+'dat/itp/'+strtrim(itp)+'/jcyl.sav'
 print,'saving the velocity field for itp=',ditp+xxx.min
 save,filename=sav_file,pror1,pror2,itp,jzt,jtt,jrt,mf_jrt,mf_jtt,mf_jzt
 

end
