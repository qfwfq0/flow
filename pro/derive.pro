function derive,path,vect1,d3
;#; this function computes the radial derivative of a input vector 
;input: vect1[*,*]: vect1[radius,modes] 
;vect1 are complex vectors


; spectrum = read_spectrum(path)
; nz = spectrum.nz
; nanf = spectrum.nanf
; mm = spectrum.mm
 mesh_size = read_mesh(path)
 dx = 1.d/double(mesh_size)
 pror1 = dindgen(mesh_size+1)/(mesh_size)
 pror2 = dindgen(mesh_size+2)/(mesh_size) - 1.d/(2.*mesh_size)

 if ( (n_elements(vect1[*,0]) gt mesh_size+2) or (n_elements(vect1[*,0]) lt mesh_size) ) then begin
  print,'errore,mesh radiale di vect1 non è giusta',n_elements(vect2[*,0])
 endif else begin
  if (n_elements(vect1[*,0]) eq mesh_size+1) then begin
;   print,'x1'
   radius = pror1
  endif else begin
;   print,'x2'
   radius = pror2
  endelse
 endelse

 array = make_array(n_elements(vect1[*,0]),n_elements(vect1[0,*]),/double,value=0.d)
 for q=0, n_elements(vect1[0,*]) -1 do begin
  dum = deriv(radius,reform(vect1[*,q]))
  array(*,q) = dum
 endfor

 return,array
end
