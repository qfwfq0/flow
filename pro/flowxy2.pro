pro flowxy, path, d3, itp, hel

RR = read_major_radius(path)
;path = './'
;d3 = 0
;hel = -10.
itp = 3999
;#; understand which kind of energy is being computed'
; restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
 mf_field='dat/ivprofiles.sav'
 dum = file_test(path+mf_field)
 print,'cacca',path+mf_field
 if (dum ne 1) then begin
  call_procedure, 'isavemodeprofiles',path,0,d3
 endif

 restore, path+mf_field
 print, 'restored!'
 print, 'helicity=', hel

 nr1 = n_elements(pror1)
 nr2 = n_elements(pror2)
 nth = nr2*2
 theta = dindgen(nth)/(nth-1) * 2. * !Dpi
 
 outfile = path+'/dat/itp/'+strtrim(string(itp,format='(i0)'))+'/flowxy.sav'
 print,'reconstruction of the flow for m=0 n=0 + m=1 n= ',hel
 ftvr = make_array(nr2,nth,/double)
 ftvt = make_array(nr2,nth,/double)
 ftvz = make_array(nr2,nth,/double)

;#; vector containing the flow in cartesian components for every time
;#; I remove the ix=-1 point of the mesh (i.e. the one at r=-0.005)
 vxy = make_array(3,nr2-1,nth,/double,value=0.)
 mm = 1
 jmn = janz(mm,hel,my,mm,nanf,nz)
 print,'jmn,',jmn
;#; axial position
 zeta = 0

 uno = make_array(n_elements(reform(real_part(pdvz[0,*,0]))),/double,value=1.)
 due = uno[1:*]#theta

;#; cycle on time to create flow
; for q = 0, n_elements(pdvz[*,0,0]) - 1 do begin
; if (q mod 10 eq 0) then print,'computing itp= ',q
;#; first: interpolation of vr on the x2 mesh
 intvr =  reform(pdvz[itp,*,*]) * dcomplex(0.d,0.d)
  for j=0, n_elements(intvr[0,*])-1 do begin
        dvr = interpol(reform(pdvr[itp,*,j]), pror1, pror2 )
        mn = mnum(j,nz,nanf,d3)
        if (mn[0] eq 0) then begin
         dvr[0] = -dvr[1]
        endif else begin
         if (mn[1] eq 1) then begin
          dvr[0] = dvr[1]
         endif else begin
          dvr[0] = -dvr[1]
         endelse
        endelse
    intvr[*,j] = dvr
  endfor
;  print,'vr interpolated'

;  for ix = 0, n_elements(pdvr[0,*,0]) - 1 do begin
;  ftvr(q,ix,*) = 2. * real_part(pdvr[q,ix,jmn]) * cos(theta) - 2. * imaginary(pdvr[q,ix,jmn]) * sin(theta)
;;#; end of cycle on ix, radial position
;  endfor  
  for ix = 0, n_elements(pdvz[0,*,0]) - 1 do begin
;#; Marco, il fattore 2. c'è già nei dati provenienti da ivprofiles.sav
   ftvr(ix,*) = real_part(intvr[ix,0])
   ftvt(ix,*) = real_part(pdvt[ix,0])
   ftvz(ix,*) = real_part(pdvz[ix,0])
;#; adding m=1 component
  for mm=1,jmn do begin
   ftvr(ix,*) =  ftvr(ix,*) + real_part(intvr[ix,jmn]) * cos(theta) -  imaginary(intvr[ix,jmn]) * sin(theta)
   ftvt(ix,*) =  ftvt(ix,*) + real_part(pdvt[itp,ix,jmn]) * cos(theta) - imaginary(pdvt[itp,ix,jmn]) * sin(theta)
   ftvz(ix,*) =  ftvz(ix,*) + real_part(pdvz[itp,ix,jmn]) * cos(theta) - imaginary(pdvz[itp,ix,jmn]) * sin(theta)
  endfor
;#; end of cycle on ix, radial position
  endfor  
  vxy[0,*,*] = ftvr[1:*,*] * cos(due) - ftvt[1:*,*] * sin(due)
  vxy[1,*,*] = ftvr[1:*,*] * sin(due) + ftvt[1:*,*] * cos(due)
  vxy[2,*,*] = ftvz[1:*,*]

;#; end of cycle on q, time
; endfor

;#; create array for partvelvec
  vvxx = make_array(n_elements(vxy[0,*,0])*n_elements(vxy[0,0,*]),/double,value=0.)
  vvyy = make_array(n_elements(vxy[0,*,0])*n_elements(vxy[0,0,*]),/double,value=0.)
  vxxx = make_array(n_elements(vxy[0,*,0])*n_elements(vxy[0,0,*]),/double,value=0.)
  vyyy = make_array(n_elements(vxy[0,*,0])*n_elements(vxy[0,0,*]),/double,value=0.)
  
;#; ordino i dati mettendo, per ogni y tutti gli x
  ny = n_elements(vxy[0,0,*])
  nx = n_elements(vxy[0,*,0])
  for q = 0, ny-1 do begin
   vvxx[q*nx:(q+1)*nx-1] = vxy[0,*,q]
   vvyy[q*nx:(q+1)*nx-1] = vxy[1,*,q]
   vxxx[q*nx:(q+1)*nx-1] = pror2[1:*]*cos(theta(q))
   vyyy[q*nx:(q+1)*nx-1] = pror2[1:*]*sin(theta(q))
  endfor
  
  save,filename=outfile,vvxx,vvyy,vxxx,vyyy
end
