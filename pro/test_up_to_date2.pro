;#; this routines check wheter my postprocessing files with the different fields are up-to-date with respect to the specyl_?_all.sav files. If not they are created/updated

;#; written by Marco Veranda, 05.09.2014

function test_up_to_date2,path,idl_file

 idl_file = path+idl_file
;#;
; str is a string with additional paramaters

;#; I count the number of specyl_?_all.sav files
 specyl_sav = file_search(path+'dat/specyl_?_all.sav',count=spcsav)

;#; I store the modification time of the higher specyl_?_all.sav file. If it is higher than the ibprofiles.sav modification time I have to rerun isavemodeprofiles.pro to obtain an up-to-date version of ibprofiles.sav
 tme_spcsav = file_info(specyl_sav(n_elements(specyl_sav)-1))
 time_spcsav = tme_spcsav.mtime 

print,'prova',idl_file
 tme_ibpr = file_info(idl_file)
 time_ibpr = tme_ibpr.mtime

 print,time_spcsav,time_ibpr
;#; I test wheter the file exists
; I test wheter ibprofiles.sav is older than any of the specyl_?_all.sav  
  if (time_spcsav gt time_ibpr) then begin
;#; need to compute it
   cmpt = 1
  endif else begin
;#; no need to compute it
   cmpt = 0
  endelse
 return, cmpt
;#;

end
