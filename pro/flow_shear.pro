pro flow_shear,path,d3,itp0,itp1,deltaitp,nth,nzz,mymin,mymax

;#; programma che calcola lo shear di velocità
;path='/ricercatori/ft/specyl/veranda/flow/archive/rfp/3d/theta16/vario_dissipazione/varioH_n7_mp4.5e-3/'
print,path

itp0 = fix(itp0,type=3)
itp1 = fix(itp1,type=3)
nitp = itp1 - itp0 + 1
nsaved = ceil(float(nitp)/deltaitp)
aitp= indgen(nsaved)*deltaitp + itp0 
print,'itp',aitp

for k=0, n_elements(aitp)-1 do begin

 iitp = aitp(k)
 anti_fft_file=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/antifft_v_'+strtrim(string(102,'(i0)'))+'x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'
 print,anti_fft_file
 anti_fft_dir=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/'
 call_procedure,'check_dir',anti_fft_dir

 ;#; check whether anti_fft_file exists
 exist = file_info(anti_fft_file)
 print,anti_fft_file
 print,exist.exists
 if (exist.exists eq 1) then begin
  restore,anti_fft_file
 endif else begin
  str = 'v'
  hel = 1.
;  nth = 128
;  nzz = 1024 
  deltaitp2 = 1
;  mymin = 0
;  mymax = 2
  ;call_procedure, 'antifft_reconst_v3', path, str, d3, hel ,nth, nzz, iitp, iitp, deltaitp2
  call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, iitp, iitp, deltaitp2, mymin, mymax
  restore,anti_fft_file
 endelse
 
 ;#; modulus of V
 vtot = sqrt(vr^2.+vt^2.+vz^2.)
 
 nr = n_elements(vtot[*,0,0])
 nt = n_elements(vtot[0,*,0])
 nz = n_elements(vtot[0,0,*])
 
 sh_tot = vtot*0.d
; grad_tot2 = vtot*0.d
; grad_tot3 = vtot*0.d
; grad_tot = vtot*0.d
 sh_r = vr*0.d
 sh_t = vt*0.d
 sh_z = vz*0.d
 for p=0,nt-1 do begin
  for q=0,nz-1 do begin
   sh_tot[*,p,q] = deriv(pror1,vtot[*,p,q])
   sh_r[*,p,q] = deriv(pror1,vr[*,p,q])
   sh_t[*,p,q] = deriv(pror2,vt[*,p,q])
   sh_z[*,p,q] = deriv(pror2,vz[*,p,q])
  endfor
 endfor
 
 nr2 = n_elements(pror2)
 shear_file=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/shear_v_'+strtrim(string(nr2,'(i0)'))+'x'+strtrim(string(nth,'(i0)'))+'_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'
 save,filename=shear_file,pror1,pror2,vecth,vecz,sh_tot,sh_r,sh_t,sh_z,iitp
; for p=1,nr-1 do begin
;  for q=0,nz-1 do begin
;   grad_tot2[p,*,q] = deriv(vecth,vtot[p,*,q])/pror1(p)
;  endfor
; endfor
; for p=1,nr-1 do begin
;  for q=0,nt-1 do begin
;   grad_tot3[p,q,*] = deriv(vecz,vtot[p,q,*])
;  endfor
; endfor
;
; ;#; gradient of velocity modulus
; grad_tot = sh_tot + grad_tot2 + grad_tot3
; ;#; relative shear
; rsh_tot = sh_tot/vtot
; rsh_r = sh_r/abs(vr)
; rsh_t = sh_t/abs(vt)
; rsh_z = sh_z/abs(vz)
 
endfor ;#; end of cicle on ITP

end   
