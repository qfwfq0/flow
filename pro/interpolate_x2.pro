function interpolate_x2, vect,path,d3

;#; interpolation of vect on the x2 mesh
;#; input: vect(*,*) : vect(radius,modes), complex
; spectrum = read_spectrum(path)
; nz = spectrum.nz
; nanf = spectrum.nanf
; mm = spectrum.mm
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
; call_procedure, 'read_spectrum',path,my,mm,nanf,nz
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
 mesh_size = read_mesh(path)
 dx = 1.d/double(mesh_size)
 pror1 = dindgen(mesh_size+1)/(mesh_size)
 pror2 = dindgen(mesh_size+2)/(mesh_size) - 1.d/(2.*mesh_size)

 ddbr =  make_array(n_elements(vect[*,0])+1,n_elements(vect[0,*]),/dcomplex)
  for j=0, n_elements(vect[0,*])-1 do begin
        dbr = interpol(reform(vect[*,j]), pror1, pror2 )
        mn = mnum(j,nz,nanf,d3)
        if (mn[0] eq 0) then begin
         dbr[0] = -dbr[1]
        endif else begin
         if (mn[1] eq 1) then begin
          dbr[0] = dbr[1]
         endif else begin
          dbr[0] = -dbr[1]
         endelse
        endelse
    ddbr[*,j] = dbr
  endfor
;print,'br interpolated'
  return,ddbr
end
