pro tau, path,d3, it_in, it_fin,tok 
;#; program computing the volume-integrated torque 
;#; computations in the FestivaldeTheorie2017 bloc notes
;#; Marco, controllato nell'aprile 2018
;#; aggiunto il torque dovuto alla sorgente di momento

; path='./'
; d3=1
; itp=200
print,'path',path
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
; call_procedure, 'read_spectrum',path,my,mm,nanf,nz,nzcon
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
 nsectionsdat = find_nsectionsdat(path)
 xxx = read_settings(path,nsectionsdat)

 RR = read_major_radius(path)
 dissipation = read_dissipation(path)
 print,'dissipation',dissipation.eta,dissipation.nu
 if (dissipation.eta*dissipation.nu le 0.d) then begin 
  print,'si è letta male la dissipazione'
  stop
 endif
; spectrum = read_spectrum(path)
 mesh_size = read_mesh(path)
 dx = 1.d/double(mesh_size)

 m0_momsour = read_momsour(path)
 print,'momsour',m0_momsour
;#; assume constant momentum source
 momsour_r = make_array(mesh_size+1,/double,value=m0_momsour(0)) 
 momsour_t = make_array(mesh_size+1,/double,value=m0_momsour(1)) 
 momsour_z = make_array(mesh_size+1,/double,value=m0_momsour(2)) 

 specyl_file = path+'dat/imf_bprofiles.sav'
 restore_bfile = path+'dat/ibprofiles.sav'
 ex_file = file_test(restore_bfile)

 if (ex_file eq 1) then begin
  part0 = test_up_to_date_files(restore_bfile,specyl_file)
  if (part0 eq 1) then begin
    call_procedure,'isavemodeprofiles',path,'0',d3
  endif
 endif

 if (ex_file ne 1) then begin
    print,'need to create the fields file ibprofiles.sav',d3
    call_procedure,'isavemodeprofiles',path,'0',d3
 endif

 part1=1 ;#; per sicurezza
 if (part1 eq 1) then begin
;#; getting the files I need
;#; I need B,v
 restore_vfile = path+'dat/ivprofiles.sav'
 if (pdbr eq !null) then begin
  restore,restore_bfile
 endif
 if (pdvr eq !null) then begin
  restore,restore_vfile
 endif

 print,'d3=',d3
; itp = 20
; for q=itp,itp do begin ;#; Marco, cycle in time

 it_in = fix(it_in)
 it_fin = fix(it_fin)
; it_fin = -1
; it_in = 250
; it_fin = 255
 if (it_fin ne -1) then begin
  ntimes = it_fin - it_in + 1
  pdtt = pdt(it_in:it_fin)
 endif else begin
  it_fin = n_elements(pdt)-1
  ntimes = it_fin-it_in+1
  pdtt = pdt(it_in:-1)
 endelse
 if (it_fin eq it_in) then begin
  ntimes = 1
  pdtt = pdt(it_fin)
 endif
 save_tau_file=path+'dat/tau_t'+strtrim(string(it_in,'(i0)'))+strtrim(string(it_fin,'(i0)'))+'.sav'
 part1 = test_up_to_date_files(save_tau_file,restore_bfile)



 tauz_jb = make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauz_vgradv =make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauz_laplv =make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauz_momsour =make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauz_tot = make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauth_jb = make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauth_vgradv =make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauth_laplv =make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauth_momsour =make_array(ntimes,mesh_size+1,/double,value=0.d)
 tauth_tot = make_array(ntimes,mesh_size+1,/double,value=0.d)

; for q=it_in,it_fin-1 do begin ;#; Marco, cycle in time
 for q=it_in,it_fin do begin ;#; Marco, cycle in time
  if (q mod 20 eq 0) then print,'time',q,pdt(q)
     ;#; per prima cosa interpolo il campo br sulla mesh x2
     ddbr = interpolate_x2(reform(pdbr[q,*,*]),path,d3) 
     ;#; per seconda cosa interpolo il campo vr sulla mesh x2
     ddvr = interpolate_x2(reform(pdvr[q,*,*]),path,d3) 

;#, first part, related to rF|theta = tau_z
;#, first term, related to r(JxB)|theta

    ;#; calcolo int_dV(brbtheta). 
     kind = 3
     parte1_1 = volume_integral(path,ddbr,reform(pdbt[q,*,*]),kind,d3,tok) 

    ;#; calcolo int_dV(rbr(dr btheta)). 
     kind = 3
     drpdbt = derive(path,reform(pdbt[q,*,*]),d3)
     rddbr = multiply_r(path,ddbr,d3)
     parte1_2 = volume_integral(path,rddbr,drpdbt,kind,d3,tok) 

    ;#; calcolo int_dV(br(dthetabr)). 
     kind = 2
     vec2 = complex( imaginary(reform(ddbr(*,*))), -real_part(reform(ddbr[*,*])) )
;     parte1_3 = volume_integral(path,ddbr,vec2,kind,d3) 
     parte1_3 = volume_integral(path,ddbr,ddbr,kind,d3,tok) 

    ;#; calcolo int_dV(bz(dthetabz)). 
     kind = 2
     vec2 = complex( imaginary(reform(pdbz(q,*,*))), -real_part(reform(pdbz[q,*,*])) )
;     parte1_4 = volume_integral(path,reform(pdbz[q,*,*]),vec2,kind,d3) 
     parte1_4 = volume_integral(path,reform(pdbz[q,*,*]),reform(pdbz[q,*,*]),kind,d3,tok) 

    ;#; calcolo int_dV(rbz(dzbtheta)). 
     kind = 1
     vec1 = multiply_r(path,reform(pdbz[q,*,*]),d3)
;     vec2 = complex( imaginary(reform(pdbt(q,*,*))), -real_part(reform(pdbt[q,*,*])) )
     vec2 = reform(pdbt[q,*,*])
     parte1_5 = volume_integral(path,vec1,vec2,kind,d3,tok) 
;  print,'fine',parte1_5

;#, second term, related to  rV.grad V|theta
    ;#; calcolo int_dV(rvr(dr vtheta)). 
     kind = 3
     drpdvt = derive(path,reform(pdvt[q,*,*]),d3)
     rddvr = multiply_r(path,ddvr,d3)
     parte2_1 = volume_integral(path,rddvr,drpdvt,kind,d3,tok) 

    ;#; calcolo int_dV(vtheta (dthetavtheta)). 
     kind = 2
     vec2 = complex( imaginary(reform(pdvt(q,*,*))), -real_part(reform(pdvt[q,*,*])) )
;     parte2_2 = volume_integral(path,reform(pdvt[q,*,*]),vec2,kind,d3) 
     parte2_2 = volume_integral(path,reform(pdvt[q,*,*]),reform(pdvt[q,*,*]),kind,d3,tok) 

    ;#; calcolo int_dV(rvz(dzvtheta)). 
     kind = 1
     vec1 = multiply_r(path,reform(pdvz[q,*,*]),d3)
;     vec2 = complex( imaginary(reform(pdvt(q,*,*))), -real_part(reform(pdvt[q,*,*])) )
     vec2 = reform(pdvt(q,*,*))
     parte2_3 = volume_integral(path,vec1,vec2,kind,d3,tok) 

    ;#; calcolo int_dV(vrvtheta). 
     kind = 3
     parte2_4 = volume_integral(path,ddvr,reform(pdvt[q,*,*]),kind,d3,tok) 

;#; third term, related to r nu laplV|theta
     
;     ;#; dr(vtheta)
;     drvtheta00 = deriv(pror2,reform(real_part(pdvt[q,*,0])))
;     dum = pror2*0.d
;     for p=1, n_elements(pror2)-1 do begin
;      temp = int_tabulated(pror2[0:p],pror2[0:p]*drvtheta00[0:p]) 
;      dum(p) = temp
;     endfor
;     parte3_1 = 4.*(!Dpi)^2.*rr(0)*dum
;
;     ;#; r d2r(vtheta)
;     d2rvtheta00 = deriv(pror2,drvtheta00)
;     dum = pror2*0.d
;     for p=1, n_elements(pror2)-1 do begin
;      temp = int_tabulated(pror2[0:p],pror2[0:p]*pror2[0:p]*d2rvtheta00[0:p]) 
;      dum(p) = temp
;     endfor
;     parte3_2 = 4.*(!Dpi)^2.*rr(0)*dum

     vt = reform(real_part(pdvt[q,*,0]))
     cc=deriv(pror2,pror2*(deriv(pror2,vt)))
     dd = pror1*0.d
     ee = pror1*0.d
     temp = 0.d
     dumm = 0.d
     for p=1, n_elements(pror1)-1 do begin
;      temp = int_tabulated(pror2[0:p],pror2[0:p]*cc[0:p]) 
;      temp = dx * pror1(p) *  0.5d * (cc(p+1)+cc(p))
;#; Marco, aprile 2018: integrale da ottenere sulla mesh X1 sfruttando quantità sulla mesh X2
      temp = temp + dx * pror2(p) * cc(p)
      dd(p) = temp
;      dumm = dx * 0.5d * (vt(p+1)+vt(p))
      dumm = dumm + dx * vt(p)
      ee(p) = dumm
     endfor
     dd = 4.d * (!Dpi)^2.d * dd * dissipation.nu 
     ee = 4.d * (!Dpi)^2.d * ee * dissipation.nu 
     

;     ;#; vtheta/r
;     dum = pror2*0.d
;     for p=1, n_elements(pror2)-1 do begin
;      temp = int_tabulated(pror2[0:p],reform(real_part(pdvt[q,0:p,0]))) 
;      dum(p) = temp
;     endfor
;     parte3_3 = 4.*(!Dpi)^2.*rr(0)*dum

     tauz_jb(q-it_in,*) = parte1_1 + parte1_2 - parte1_3 - parte1_4 + parte1_5
     tauz_vgradv(q-it_in,*) = parte2_1 + parte2_2 + parte2_3 + parte2_4
;     tauz_laplv(q,*) = (parte3_1 + parte3_2 - parte3_3) * dissipation.nu
     tauz_laplv(q-it_in,*) = (dd - ee) 
     tauz_momsour(q-it_in,*) = pror1 * momsour_t

     tauz_tot(q-it_in,*) = tauz_jb(q-it_in,*) + tauz_vgradv(q-it_in,*) + tauz_laplv(q-it_in,*) + tauz_momsour(q-it_in,*)

;#, second part, related to rF|z = -tau_theta
;#, first term, related to r(JxB)|z

    ;#; calcolo int_dV(btheta(dthetabz)). 
     kind = 2
;     vec2 = complex( imaginary(reform(pdbz(q,*,*))), -real_part(reform(pdbz[q,*,*])) )
;     parte4_1 = volume_integral(path,reform(pdbt[q,*,*]),vec2,kind,d3) 
     parte4_1 = volume_integral(path,reform(pdbt[q,*,*]),reform(pdbz[q,*,*]),kind,d3,tok) 

    ;#; calcolo int_dV(rbtheta(dzbtheta)). 
     kind = 1
     vec1 = multiply_r(path,reform(pdbt[q,*,*]),d3)
;     vec2 = complex( imaginary(reform(pdbt(q,*,*))), -real_part(reform(pdbt[q,*,*])) )
     vec2 = reform(pdbt[q,*,*])
     parte4_2 = volume_integral(path,vec1,vec2,kind,d3,tok) 

    ;#; calcolo int_dV(rbr(dzbr)). 
     kind = 1
     vec1 = multiply_r(path,reform(ddbr[*,*]),d3)
;     vec2 = complex( imaginary(reform(ddbr(*,*))), -real_part(reform(ddbr[*,*])) )
     vec2 = reform(ddbr(*,*))
     parte4_3 = volume_integral(path,vec1,vec2,kind,d3,tok) 

    ;#; calcolo int_dV(rbr(drbz)). 
     kind = 3
     vec1 = multiply_r(path,reform(ddbr[*,*]),d3)
     drpdbz = derive(path,reform(pdbz[q,*,*]),d3)
     parte4_4 = volume_integral(path,vec1,drpdbz,kind,d3,tok) 

;#, second term, related to  rV.grad V|z
    ;#; calcolo int_dV(rvr(dr vz)). 
     kind = 3
     drpdvz = derive(path,reform(pdvz[q,*,*]),d3)
     rddvr = multiply_r(path,ddvr,d3)
     parte5_1 = volume_integral(path,rddvr,drpdvz,kind,d3,tok) 

    ;#; calcolo int_dV(vtheta (dthetavz)). 
     kind = 2
;     vec2 = complex( imaginary(reform(pdvz(q,*,*))), -real_part(reform(pdvz[q,*,*])) )
;     parte5_2 = volume_integral(path,reform(pdvt[q,*,*]),vec2,kind,d3) 
     parte5_2 = volume_integral(path,reform(pdvt[q,*,*]),reform(pdvz[q,*,*]),kind,d3,tok) 

    ;#; calcolo int_dV(rvz(dzvz)). 
     kind = 1
     vec1 = multiply_r(path,reform(pdvz[q,*,*]),d3)
;     vec2 = complex( imaginary(reform(pdvz(q,*,*))), -real_part(reform(pdvz[q,*,*])) )
;     parte5_3 = volume_integral(path,vec1,vec2,kind,d3) 
     parte5_3 = volume_integral(path,vec1,reform(pdvz[q,*,*]),kind,d3,tok) 

;#; third term, related to r nu laplV|z

;     ;#; dr(vz)
;     drvz00 = deriv(pror2,reform(real_part(pdvz[q,*,0])))
;     dum = pror2*0.d
;     for p=1, n_elements(pror2)-1 do begin
;      temp = int_tabulated(pror2[0:p],pror2[0:p]*drvz00[0:p]) 
;      dum(p) = temp
;     endfor
;     parte6_1 = 4.*(!Dpi)^2.*rr(0)*dum
;
;     ;#; r d2r(vz)
;     d2rvz00 = deriv(pror2,drvz00)
;     dum = pror2*0.d
;     for p=1, n_elements(pror2)-1 do begin
;      temp = int_tabulated(pror2[0:p],pror2[0:p]*pror2[0:p]*d2rvz00[0:p]) 
;      dum(p) = temp
;     endfor
;     parte6_2 = 4.*(!Dpi)^2.*rr(0)*dum
     vz = reform(real_part(pdvz[q,*,0]))
     cc=deriv(pror2,pror2*(deriv(pror2,vz)))
     dd = pror1*0.d
;     ee = pror1*0.d
     temp = 0.d
     for p=1, n_elements(pror1)-1 do begin
;      temp = int_tabulated(pror2[0:p],pror2[0:p]*cc[0:p]) 
;      temp = dx * pror1(p) *  0.5d * (cc(p+1)+cc(p))
;#; Marco, aprile 2018: integrale da ottenere sulla mesh X1 sfruttando quantità sulla mesh X2
      temp = temp + dx * pror2(p) *  cc(p) ;0.5d * (cc(p+1)+cc(p))
      dd(p) = temp
     endfor
     dd = 4.d*(!Dpi)^2.d * dd * dissipation.nu 

;#; i segni meno davanti sono dovuti al fatto che tau_theta = r x F = -rF|z
     tauth_jb(q-it_in,*) = -(parte4_1 - parte4_2 - parte4_3 + parte4_4)
     tauth_vgradv(q-it_in,*) = -(parte5_1 + parte5_2 + parte5_3)
;     tauth_laplv(q,*) = -(parte6_1 + parte6_2 ) * dissipation.nu
     tauth_laplv(q-it_in,*) = - dd
     tauth_momsour(q-it_in,*) = - pror1 * momsour_z

     tauth_tot(q-it_in,*) = tauth_jb(q-it_in,*) + tauth_vgradv(q-it_in,*) + tauth_laplv(q-it_in,*) + tauth_momsour(q-it_in,*)
 endfor ;#; Marco, end of cycle in time
    
 save,pdtt,pror1,pror2,tauz_jb,tauz_vgradv,tauz_laplv,tauz_momsour,tauz_tot,tauth_jb,tauth_vgradv,tauth_laplv,tauth_momsour,tauth_tot,filename=save_tau_file

 endif else begin ;#; end of part that computed volume integrated tau
  ;#; if save_tau_file already exists, simply restore it
  print,'I restore the already existing tau file'
  restore,save_tau_file
 endelse

 save_angularmom_file=path+'dat/angular_momentum_t'+strtrim(string(it_in,'(i0)'))+strtrim(string(it_fin,'(i0)'))+'.sav'
 part1 = test_up_to_date_files(save_angularmom_file,save_tau_file)
 if (part1 eq 1) then begin
  print, 'I compute the angular momentum'
;
  ang_z = tauz_tot * 0.d
  ang_th = tauth_tot * 0.d
  for p = 0, n_elements(tauz_tot[0,*])-1 do begin ;#; cycle on radius
   for q = 1, n_elements(tauz_tot[*,0])-1 do begin ;#; cycle on time
;    print,'q',q
    ang_z(q,p) = ang_z(q-1,p) + (pdt(q)-pdt(q-1)) * tauz_tot(q,p)
    ang_th(q,p) = ang_th(q-1,p) + (pdt(q)-pdt(q-1)) * tauth_tot(q,p)
   endfor
  endfor
  save,pdtt,pror1,pror2,tauz_tot,tauth_tot,ang_z,ang_th,filename=save_angularmom_file
;  
 endif else begin
  print,'ho già calcolato tutto!!! tau, momento angolare'
 endelse
end

