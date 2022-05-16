;#; written by Marco Veranda, 24.11.2015
;#; this routine checks wheter file1 is older than file2.
;#; if (file1 older than file2) need to recompute file 1, so cmpt=1

function test_up_to_date_files,idl_file1,idl_file2

 idl_file1 = idl_file1
 idl_file2 = idl_file2
; print,'a',idl_file1
; print,'b',idl_file2
;#;

;#; I store the modification time of the higher specyl_?_all.sav file. If it is higher than the ibprofiles.sav modification time I have to rerun isavemodeprofiles.pro to obtain an up-to-date version of ibprofiles.sav
 tme_sav1 = file_info(idl_file1)
 time_sav1 = tme_sav1.mtime 

 tme_sav2 = file_info(idl_file2)
 time_sav2 = tme_sav2.mtime

 print,time_sav1,time_sav2
;#; I test wheter the file exists
; I test wheter ibprofiles.sav is older than any of the specyl_?_all.sav  
  if (time_sav1 lt time_sav2) then begin
;#; need to compute it
   cmpt = 1
  endif else begin
;#; no need to compute it
   cmpt = 0
  endelse
 return, cmpt
;#;

end
