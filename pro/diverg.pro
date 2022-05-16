function diverg,path,vr,vt,vz,j,d3
;#; program computing divergence of a field in cylindrical geometry

;path='./'
;j=73
;d3=1
;restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
restore,path+'dat/ibprofiles.sav'
RR = read_major_radius(path)

dx = pror1(5)-pror1(4)

;
;vr = reform(pdbr[100,*,j])
;vt = reform(pdbt[100,*,j])
;vz = reform(pdbz[100,*,j])

;#; 1/r dr(r*vr)
dvr = vt*dcomplex(0.d,0.d)
for i = 1, n_elements(dvr)-2 do begin
 dvr(i) = (pror1(i)*vr(i) - pror1(i-1)*vr(i-1)) / dx
endfor
 mn = mnum(j,nz,nanf,d3)
;#; regularity at edge
  dvr(-1) = dvr(n_elements(dvr)-2)
;#; regularity at core
 if (mn[0] eq 0) then begin
  dvr(0) = + dvr(1)
  vt(0) = - vt(1)
 endif else begin
  if (mn[0] eq 1) then begin
   dvr(0) = - dvr(1)
   vt(0) = + vt(1)
  endif else begin
   dvr(0) = + dvr(1)
  vt(0) = - vt(1)
  endelse
 endelse

 ddvr = dvr / pror2

 ddvt = mn[0] * vt / pror2
 ddvt(0) = + ddvt(1)
 ddvz = mn[1] / RR[0] * vz

 rdiv = real_part(ddvr) - imaginary(ddvt) - imaginary(ddvz)
 idiv = imaginary(ddvr) + real_part(ddvt) + real_part(ddvz)
 div = dcomplex(rdiv,idiv)
 return,div
end
