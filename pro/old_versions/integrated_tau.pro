pro integrated_tau, path,d3
;#; program computing the momentum of the JxB force in cylindrical geometry, z component
;#; Marco, added also
;#; program computing the JxB force in cylindrical geometry, z component
;#; Marco, 14/5/2015 added momentum of the JxB force, th component

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

; itp = xxx.max
; help,xxx,/str
 
; if (itp lt xxx.min) then begin
;  print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
;  stop
; endif
 
 RR = read_major_radius(path)
;#; getting the files I need
;#; I need B
 restore_ifile = path+'dat/ibprofiles.sav'
 ex_file = file_test(restore_ifile)
 if (ex_file ne 1) then begin
    print,'need to create the fields file ibprofiles.sav',d3
    call_procedure,'isavemodeprofiles',path,'0',d3
 endif
 restore,restore_ifile
 ddr = pror2(5) - pror2(4)
;#; First of all I need B^2
; restore_ifile = path+'dat/ibsq_profiles.sav'
; ex_file = file_test(restore_ifile)
; if (ex_file ne 1) then begin
;  print,'need to create the fields file ibsq_profiles.sav'
;  call_procedure,'squared',path,'b',d3
; endif
; restore,restore_ifile

; ditp=fix(itp-xxx.min) 
; square = reform(square(ditp,*,*))
; pdbr = reform(pdbr[ditp,*,*])
; pdbt = reform(pdbt[ditp,*,*])
; pdbz = reform(pdbz[ditp,*,*])
; tau_th = pdbz * dcomplex(0.d,0.d)
; tau_z = pdbz * dcomplex(0.d,0.d)
 brbth = make_array(n_elements(pdbt[0,*,0]),/double,value=0.) 
 brbz = make_array(n_elements(pdbt[0,*,0]),/double,value=0.) 
 momth_1 = make_array(n_elements(pdbt[0,*,0]),/double,value=0.) 
 momth_2 = make_array(n_elements(pdbt[0,*,0]),/double,value=0.) 
 momth_dum = make_array(n_elements(pdbt[0,*,0]),/double,value=0.) 
 int_tau_z = make_array(n_elements(pdbt[*,0,0]),n_elements(pdbt[0,*,0]),/double,value=0.) 
 int_tau_th = make_array(n_elements(pdbt[*,0,0]),n_elements(pdbt[0,*,0]),/double,value=0.) 
 int_jb_z = make_array(n_elements(pdbt[*,0,0]),n_elements(pdbt[0,*,0]),/double,value=0.) 

;#; component z
;#; tauz = -0.5 * dth(B^2) + div(rr*Bth*B)
;#; component th
;#; tauth = 0.5 * rr * dth(B^2) - rr * div(Bth*B)
;#; volume integrated component z
;#; int_tauz / (4*pi^2*R0)= 4. * rr^2 * (Re(br)Re(bth) + Im(br)Im(bth))
 dthb2 = pdbt * dcomplex(0.d,0.d)
 dzb2  = pdbt * dcomplex(0.d,0.d)

;#; first: interpolation of br on the x2 mesh
 ddbr =  pdbt * dcomplex(0.d,0.d)
 help,pdbt
 print,nanf
  for q=0, n_elements(pdbt[*,0,0])-1 do begin
  for j=0, n_elements(pdbt[0,0,*])-1 do begin
        dbr = interpol(reform(pdbr[q,*,j]), pror1, pror2 )
        mn = mnum(j,nz,nanf,d3)
        if (mn[0] eq 0) then begin
         dbr[0] = -dbr[1]
        endif else begin
         if (mn[0] eq 1) then begin
          dbr[0] = dbr[1]
         endif else begin
          dbr[0] = -dbr[1]
         endelse
        endelse
    ddbr[q,*,j] = dbr
  endfor
  endfor
print,'br interpolated'

  for q=0, n_elements(pdbt[*,0,0])-1 do begin
  if (q mod 50 eq 0) then print,'computed itp=',q
  brbth = brbth * 0.d
  brbz = brbz * 0.d
  momth_2 = momth_2 * 0.d
  for j=0, n_elements(pdbt[0,0,*])-1 do begin
;#; seventh: computation of eq16 of Pustovitov NF47, 1583 (2007)
;   rebrrebth = nonlinearconvo(path,reform(real_part(ddbr)),reform(real_part(pdbt)),j,d3)
;   imbrimbth = nonlinearconvo(path,reform(imaginary(ddbr)),reform(imaginary(pdbt)),j,d3)

;#; I do not consider j=0 because, up to now, all the computations involve br (that is zero in equilibrium)
    if (j gt 0) then begin
     if (d3 eq 0) then begin
      brbth = brbth + real_part(ddbr[q,*,j])*real_part(pdbt[q,*,j]) + imaginary(ddbr[q,*,j])*imaginary(pdbt[q,*,j])
;#; this quantity is needed also for the theta component of the momentum
      brbz = brbz + real_part(ddbr[q,*,j])*real_part(pdbz[q,*,j]) + imaginary(ddbr[q,*,j])*imaginary(pdbz[q,*,j])
     endif else begin
      print,'3d case not implemented yet'
     endelse
;#; for the theta component I need an integral in the radial direction
;#; attention to the first cell of the integral
    momth_dum = momth_dum * 0.d
    for i=1, n_elements(pdbt[0,*,0])-1 do begin
     if (i eq 1) then begin
      momth_dum[i] = ddr / 2.d * ( pror2(i) * (real_part(ddbr[q,i,j])*real_part(pdbz[q,i,j]) + imaginary(ddbr[q,i,j])*imaginary(pdbz[q,i,j])) )
     endif else begin
      momth_dum[i] = momth_dum[i-1] + ddr * ( pror2(i) * (real_part(ddbr[q,i,j])*real_part(pdbz[q,i,j]) + imaginary(ddbr[q,i,j])*imaginary(pdbz[q,i,j])) )
     endelse
    endfor
    momth_2 = momth_2 + momth_dum

    endif

   int_tau_th[q,*] = - 4.d *  ( pror2^2.d * brbz - momth_2)
   int_tau_z[q,*] = 4.d * pror2^2.d * brbth
   int_jb_z[q,*] = 4.d * pror2 * brbz
 endfor
 endfor


;#; saving 
; sav_file = path+'dat/itp/'+strtrim(itp)+'/tau.sav'
 sav_file = path+'dat/integrated_tau.sav'
 print,'saving'
 save, filename=sav_file, pdt,pror2,int_tau_z,int_jb_z,int_tau_th
end
