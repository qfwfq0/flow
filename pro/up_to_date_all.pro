;#; this routines check wheter my postprocessing files specyl_?_all.sav are up-to-date with respect to the specyl_?_all.dat files. If not they are created/updated

;#; written by Marco Veranda, 12.02.2015

function up_to_date_all,path,nsec

;#; I count the number of specyl_?_all.sav files
 specyl_sav = file_search(path+'dat/specyl_'+strtrim(string(nsec,format='(i0)'))+'_all.sav',count=spcsav)
 specyl_dat = file_search(path+'dat/specyl_'+strtrim(string(nsec,format='(i0)'))+'_all.dat',count=spcdat)

  tme_spcsav = file_info(specyl_sav)
  time_spcsav = tme_spcsav.mtime 
  tme_spcdat = file_info(specyl_dat)
  time_spcdat = tme_spcdat.mtime

  if (time_spcsav lt time_spcdat) then begin
   compute_it = 1
  endif else begin
   compute_it = 0
  endelse
;;#; I test wheter the file exists
; exist_file = file_test(idl_file)
; if (exist_file ne 1) then begin
;   print,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
;   print,'Creating  '+call_proc
;   print,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
;  call_procedure,call_proc,path
; endif else begin
;; I test wheter ibprofiles.sav is older than any of the specyl_?_all.sav  
;  if (time_spcsav gt time_ibpr) then begin
;   print,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
;   print,'Updating  '+call_proc
;   print,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
;   call_procedure,call_proc,path
;  endif else begin
;  endelse
; endelse

 return,compute_it
;#;

end
