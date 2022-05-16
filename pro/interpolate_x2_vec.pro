function interpolate_x2_vec, vect,path,d3

;#; interpolation of vect on the x2 mesh
;#; input: vect(*,*,*) : vect(radius,theta,zeta), complex
; spectrum = read_spectrum(path)
; nz = spectrum.nz
; nanf = spectrum.nanf
; mm = spectrum.mm
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
 nzcon=nz
 nzcon[0]=nz[0]+1
endelse
 mesh_size = read_mesh(path)
 dx = 1.d/double(mesh_size)
 pror1 = dindgen(mesh_size+1)/(mesh_size)
 pror2 = dindgen(mesh_size+2)/(mesh_size) - 1.d/(2.*mesh_size)

 ddbr =  make_array(n_elements(vect[*,0,0])+1,n_elements(vect[0,*,0]),n_elements(vect[0,0,*]),/dcomplex)
  for j=0, n_elements(vect[0,*,0])-1 do begin
  for p=0, n_elements(vect[0,0,*])-1 do begin
    dbr = interpol(reform(vect[*,j,p]), pror1, pror2 )
    ddbr[*,j,p] = dbr
  endfor
  endfor
;print,'br interpolated'
  return,ddbr
end
