pro idl_saveb, path, itp

 nsectionsdat = find_nsectionsdat(path)
 xxx = read_settings(path,nsectionsdat)

 print,itp
 help,xxx,/str
 idl_file=path+'dat/ibprofiles.sav'
 print,'restoring...',idl_file
 restore,idl_file
 
 if (itp lt xxx.min) then begin
  print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
  stop
 endif
 
 ditp=fix(itp-xxx.min) 
 bzt=make_array(n_elements(pdbz[0,*,0]),n_elements(pdbz[0,0,*]),/dcomplex,value=0.)
 btt=make_array(n_elements(pdbt[0,*,0]),n_elements(pdbt[0,0,*]),/dcomplex,value=0.)
 brt=make_array(n_elements(pdbr[0,*,0]),n_elements(pdbr[0,0,*]),/dcomplex,value=0.)
 
 bzt[*,*]=reform(pdbz[ditp,*,*]) 
 btt[*,*]=reform(pdbt[ditp,*,*]) 
 brt[*,*]=reform(pdbr[ditp,*,*]) 
 
;#; write the field using the modulus-phase convention
 mf_bzt=make_array(2,n_elements(pdbz[0,*,0]),n_elements(pdbz[0,0,*]),/double,value=0.)
 mf_btt=make_array(2,n_elements(pdbt[0,*,0]),n_elements(pdbt[0,0,*]),/double,value=0.)
 mf_brt=make_array(2,n_elements(pdbr[0,*,0]),n_elements(pdbr[0,0,*]),/double,value=0.)
 mf_bzt[0,*,*] = sqrt( real_part(bzt[*,*])^2.d + imaginary(bzt[*,*])^2.d )
 mf_btt[0,*,*] = sqrt( real_part(btt[*,*])^2.d + imaginary(btt[*,*])^2.d )
 mf_brt[0,*,*] = sqrt( real_part(brt[*,*])^2.d + imaginary(brt[*,*])^2.d )
 mf_bzt[1,*,*] = atan( imaginary(bzt[*,*]),real_part(bzt[*,*]) )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_bzt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_bzt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_bzt[1,aaa[0,k],aaa[1,k]] = mf_bzt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

 mf_btt[1,*,*] = atan( imaginary(btt[*,*]) , real_part(btt[*,*]) )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_btt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_btt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_btt[1,aaa[0,k],aaa[1,k]] = mf_btt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

 mf_brt[1,*,*] = atan( imaginary(brt[*,*]), real_part(brt[*,*]) )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_brt[1,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_brt[1,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_brt[1,aaa[0,k],aaa[1,k]] = mf_brt[1,aaa[0,k],aaa[1,k]] + 2.*!dpi
  endfor
 endif

;#; saving 
;#; check directory existence
 dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
 call_procedure,'check_dir',dir 
  
 sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/bcyl.sav'
 print,'saving the velocity field for itp=',ditp+xxx.min
 save,filename=sav_file,pror1,pror2,itp,bzt,btt,brt,mf_brt,mf_btt,mf_bzt
 

end
