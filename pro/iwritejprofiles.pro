pro iwritejprofiles,path,mf,d3,itp0

 save_file = path+'dat/ijprofiles.sav'
 result = file_test(save_file)
 mf = 0
 
 if (result ne 1) then begin
  call_procedure, 'isavejprofiles', path, mf, d3
 endif else begin
  restore,save_file
 endelse
 
 print, itp0
 current_r = reform(pjr(itp0,*,*))
 current_t = reform(pjt(itp0,*,*))
 current_z = reform(pjz(itp0,*,*))

 call_procedure,'write_to_fortran',path,pror1,pror2,current_r,current_t,current_z,itp0
end
