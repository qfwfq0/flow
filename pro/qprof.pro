pro qprof, path

   RR = read_major_radius(path)
 mf_field='dat/imf_bprofiles.sav'
 dum = file_test(path+mf_field)
 print,'testing wheter imf_bprofiles.sav exists',path+mf_field
 if (dum ne 1) then begin
  nsectionssav = find_nsectionssav(path)
  mf = 1
  call_procedure, 'ireadsc2',path,nsectionssav,mf,d3
 endif

 restore, path+mf_field
 print, 'restored!'

 outfile = path+'/dat/qprof.sav'
;#; if cycle on the file time_chihel.sav existence
 dum = file_test(outfile)
 tme_sav = file_info(outfile)
 time_sav = tme_sav.mtime 
 tme_profsav = file_info(path+mf_field)
 time_profsav = tme_profsav.mtime 
 if (dum ne 1 or time_sav lt time_profsav) then begin
  uno = make_array(n_elements(reform(bz[*,0,0,0])),/double,value=1.)
  r2r2 = uno#pror2
  qprof=bz[*,*,0,0]*cos(bz[*,*,0,1])*r2r2/(4.d0*bt[*,*,0,0]*cos(bt[*,*,0,1]))
;#; average of qprof between radius 0 and radius 2
  q02 = make_array(n_elements(qprof[*,0]),/double,value=0.) 
  q05 = make_array(n_elements(qprof[*,0]),/double,value=0.) 
  q010 = make_array(n_elements(qprof[*,0]),/double,value=0.) 
  q95 = make_array(n_elements(qprof[*,0]),/double,value=0.) 
  dqdr = make_array(n_elements(qprof[*,0]),n_elements(qprof[0,*]),/double,value=0.) 
  shear_q = make_array(n_elements(qprof[*,0]),n_elements(qprof[0,*]),/double,value=0.) 
  for p=0,n_elements(qprof[*,0])-1 do begin
   aaa=moment(qprof[p,0:2])
   q02[p] = aaa[0]
   aaa=moment(qprof[p,0:5])
   q05[p] = aaa[0]
   aaa=moment(qprof[p,0:10])
   q010[p] = aaa[0]
   q95[p] = reform(qprof[p,-5])
   dqdr[p,*] =  deriv(pror2,qprof[p,*])
   shear_q[p,*] =  deriv(pror2,qprof[p,*]) * (pror2/qprof[p,*])
  endfor

;#; save
 save,filename=outfile,qprof,tt,pror2,q02,q05,q010,q95,dqdr,shear_q
 endif else begin
;#; restore already computed
 restore,filename=outfile
 endelse

end

