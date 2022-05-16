;#; written by Marco Veranda, 08.09.2016
;#; this routine removes useless *sav files if the simulation is ended

pro clean_specyl,path

;#; I count the number of specyl_?_all.sav files
 specyl_sav = file_search(path+'dat/specyl_?_all.sav',count=spcsav)
;#; I store the modification time of the higher specyl_?_all.sav file.
 tme_spc = file_info(specyl_sav(spcsav-1))
 time_spc = tme_spc.mtime 
 
;#; find imf_bprofiles_?.sav files
 imf_b = file_search(path+'dat/imf_bprofiles_?.sav',count=imf_bc)
;#; I store the modification time of the higher imf_b_profiles_?.sav file.
 tme_imf_b = file_info(imf_b(imf_bc-1))
 time_imf_b = tme_imf_b.mtime 

 print,'     '
 if (time_imf_b gt time_spc) then begin
  print,'delete the imf_bprofiles_?.sav files'
  for i=0, n_elements(imf_b)-1 do begin
   print,imf_b(i)
;   spawn,'rm -f '+imf_b(i)
  endfor
 endif

 print,'     '
;#; find imf_bprofiles_?.sav files
 imf_v = file_search(path+'dat/imf_vprofiles_?.sav',count=imf_vc)
;#; I store the modification time of the higher imf_v_profiles_?.sav file.
 tme_imf_v = file_info(imf_v(imf_vc-1))
 time_imf_v = tme_imf_v.mtime 

 if (time_imf_v gt time_spc) then begin
  print,'delete the imf_vprofiles_?.sav files'
  for i=0, n_elements(imf_v)-1 do begin
   print,imf_v(i)
;   spawn,'rm -f '+imf_b(i)
  endfor
 endif

end
