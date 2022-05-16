pro idl_savev, path, itp

 nsectionsdat = find_nsectionsdat(path)
 xxx = read_settings(path,nsectionsdat)

 print,itp
 help,xxx,/str
 idl_file=path+'dat/ivprofiles.sav'
 print,'restoring...',idl_file
 restore,idl_file
 
 if (itp lt xxx.min) then begin
  print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
  stop
 endif
 
 ditp=fix(itp-xxx.min) 
 vzt=make_array(n_elements(pdvz[0,*,0]),n_elements(pdvz[0,0,*]),/dcomplex,value=0.)
 vtt=make_array(n_elements(pdvt[0,*,0]),n_elements(pdvt[0,0,*]),/dcomplex,value=0.)
 vrt=make_array(n_elements(pdvr[0,*,0]),n_elements(pdvr[0,0,*]),/dcomplex,value=0.)
 
 vzt[0,*,*]=reform(pdvz[ditp,*,*]) 
 vtt[0,*,*]=reform(pdvt[ditp,*,*]) 
 vrt[0,*,*]=reform(pdvr[ditp,*,*]) 
 
;#; write the field using the modulus-phase convention
 mf_vzt=make_array(2,n_elements(pdvz[0,*,0]),n_elements(pdvz[0,0,*]),/double,value=0.)
 mf_vtt=make_array(2,n_elements(pdvt[0,*,0]),n_elements(pdvt[0,0,*]),/double,value=0.)
 mf_vrt=make_array(2,n_elements(pdvr[0,*,0]),n_elements(pdvr[0,0,*]),/double,value=0.)
 mf_vzt[0,*,*] = sqrt( real_part(vzt[*,*]^2.d) + imaginary(vzt[*,*])^2.d )
 mf_vtt[0,*,*] = sqrt( real_part(vtt[*,*]^2.d) + imaginary(vtt[*,*])^2.d )
 mf_vrt[0,*,*] = sqrt( real_part(vrt[*,*]^2.d) + imaginary(vrt[*,*])^2.d )
 mf_vzt[1,*,*] = atan( imaginary(vzt[*,*]), real_part(vzt[0,*,*]) )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_vzt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_vzt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_vzt[1,aaa[0,k],aaa[1,k]] = mf_vzt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

 mf_vtt[1,*,*] = atan( imaginary(vtt[*,*]), real_part(vtt[*,*]) )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_vtt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_vtt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_vtt[1,aaa[0,k],aaa[1,k]] = mf_vtt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

 mf_vrt[1,*,*] = atan( imaginary(vrt[*,*]) ,real_part(vrt[*,*]) )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_vrt[1,*,*]) lt 0.)
  if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_vrt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_vrt[1,aaa[0,k],aaa[1,k]] = mf_vrt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

;#; save
 sav_file = path+'dat/itp/'+strtrim(itp)+'/vcyl.sav'
 print,'saving the velocity field for itp=',ditp+xxx.min
 save,filename=sav_file,pror1,pror2,itp,vzt,vtt,vrt,mf_vrt,mf_vtt,mf_vzt
end
