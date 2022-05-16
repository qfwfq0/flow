pro momentum, path,d3
;#; program computing the volume-integrated momentum in z direction
;#; computations in the FestivaldeTheorie2017 bloc notes

; path='./'
; d3=1
; itp=200
; restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse

 nsectionsdat = find_nsectionsdat(path)
 xxx = read_settings(path,nsectionsdat)

; if (itp lt xxx.min) then begin
;  print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
;  stop
; endif
; 
 RR = read_major_radius(path)
 dissipation = read_dissipation(path)
; print,'dissipation',dissipation
 spectrum = read_spectrum(path)
 mesh_size = read_mesh(path)
; print,'spectrum',spectrum

;#; getting the files I need
;#; I need B,v
 restore_ifile = path+'dat/ibprofiles.sav'
 restore_vfile = path+'dat/ivprofiles.sav'
 ex_file = file_test(restore_ifile)
 if (ex_file ne 1) then begin
    print,'need to create the fields file ibprofiles.sav',d3
    call_procedure,'isavemodeprofiles',path,'0',d3
 endif
 if (pdbr eq !null) then begin
  restore,restore_ifile
 endif
 if (pdvr eq !null) then begin
  restore,restore_vfile
 endif

 print,'d3=',d3
; itp = 20
; for q=itp,itp do begin ;#; Marco, cycle in time

 it_in = 0
 it_fin = -1
 if (it_fin ne -1) then begin
  ntimes = it_fin - it_in + 1
 endif else begin
  it_fin = n_elements(pdt)
  ntimes = it_fin - it_in + 1
 endelse
 pz_tot = make_array(ntimes,mesh_size+2,/double,value=0.d)
 pth_tot = make_array(ntimes,mesh_size+2,/double,value=0.d)

 for q=it_in,it_fin-1 do begin ;#; Marco, cycle in time

  if (q mod 20 eq 0) then print,'time',q,pdt(q)
;#, first term, related to JxB|z

    ;#; intanto calcolo int_dV(dz(b^2)). Sui motivi per cui passo alla subroutine proprio quel valore di vec2, gaurdare sul bloc notes AIX2017 (è dovuto alla derivata in zeta, che mischia parti reali e immaginarie del vettore derivato). Analogo per i casi successivi.
     kind = 1
     vec1 = reform(pdbz(q,*,*))
     vec2 = complex( imaginary(reform(pdbz(q,*,*))), -real_part(reform(pdbz[q,*,*])) )
     parte1 = volume_integral(path,vec1,vec2,kind,d3) 
     ;#; per prima cosa interpolo il campo br sulla mesh x2
     ddbr = interpolate_x2(reform(pdbr[q,*,*]),path,d3) 
     vec1 = reform(ddbr(*,*))
     vec2 = complex( imaginary(reform(ddbr[*,*])), -real_part(reform(ddbr[*,*])) )
     parte2 = volume_integral(path,vec1,vec2,kind,d3) 
     vec1 = reform(pdbt(q,*,*))
     vec2 = complex( imaginary(reform(pdbt(q,*,*))), -real_part(reform(pdbt[q,*,*])) )
     parte3 = volume_integral(path,vec1,vec2,kind,d3) 
     first = 2.d * (parte1 + parte2 + parte3)

    ;#; continuo con int_dS3(br*bz)) (integrale su una superficie a guscio cilindrico)
     parte4 = surface_integral(path,ddbr,reform(pdbz[q,*,*]),kind,d3) 

;#, second term, related to V.grad V|z
     ;#; vr*dr(vz)
     ;#; devo creare una routine che fa la derivata radiale
     kind = 3
     ;#; per prima cosa interpolo il campo vr sulla mesh x2
     ddvr = interpolate_x2(reform(pdvr[q,*,*]),path,d3) 
     dvr = derive(path,ddvr,d3)
     parte5 = volume_integral(path,ddvr,dvr,kind,d3) 
     ;#; vt/r * d(theta)vz
     kind = 4
     vec1 = reform(pdvt(q,*,*))
     vec2 = complex( imaginary(reform(pdvz(q,*,*))), -real_part(reform(pdvz[q,*,*])) )
     parte6 = volume_integral(path,vec1,vec2,kind,d3) 

     kind = 1
     vec1 = reform(pdvz(q,*,*))
     vec2 = complex( imaginary(reform(pdvz(q,*,*))), -real_part(reform(pdvz[q,*,*])) )
     parte7 = volume_integral(path,vec1,vec2,kind,d3) 

;#; third term, related to nu.laplV|z
    
     drvz00 = deriv(pror2,reform(real_part(pdvz[q,*,0])))
     parte8 = 4.*(!Dpi)^2.*rr(0)*dissipation.nu*pror2*drvz00

;;Marco, attenzione ai segni delle varie quantità
     pz_jb(q,*) = - (parte1 + parte2 +parte3) + parte4 
     pz_vgradv(q,*) = parte4 - (parte5 + parte6 + parte7)
     pz_laplv(q,*) = parte8
     pz_tot(q,*) = - (parte1 + parte2 +parte3) + parte4 - ( parte5 + parte6 +parte7 ) + parte8

;###############;
;#; inizia il calcolo del momento in direzione theta
;###############;

;#, first term, related to JxB|theta

    ;#; intanto calcolo -1/(2*r)*int_dV(dtheta(b^2)). Sui motivi per cui passo alla subroutine proprio quel valore di vec2, gaurdare sul bloc notes AIX2017 (è dovuto alla derivata in theta, che mischia parti reali e immaginarie del vettore derivato). Analogo per i casi successivi.
     ;#; 1/r * br*d(theta)br
;     kind = 4
;     vec2 = complex( imaginary(reform(ddbr(*,*))), -real_part(reform(ddbr[*,*])) )
;     theta1 = volume_integral(path,ddbr,vec2,kind,d3) 
;
;     ;#; 1/r * bt*d(theta)bt
;     kind = 4
;     vec2 = complex( imaginary(reform(pdbt(q,*,*))), -real_part(reform(pdbt[q,*,*])) )
;     theta2 = volume_integral(path,reform(pdbt(q,*,*)),vec2,kind,d3) 
;
;     ;#; 1/r * bz*d(theta)bz
;     kind = 4
;     vec2 = complex( imaginary(reform(pdbz(q,*,*))), -real_part(reform(pdbz[q,*,*])) )
;     theta3 = volume_integral(path,reform(pdbz(q,*,*)),vec2,kind,d3) 

 endfor ;#; Marco, end of cycle in time
    
 save_file=path+'dat/imomentum.sav'
 save,pdt,pror1,pror2,pz_tot,filename=save_file
stop
end
