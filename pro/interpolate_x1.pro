function interpolate_x1, vect

;#; interpolation of vect on the x1 mesh
;#; input: vect(*,*,*) : vect(radius,theta,zeta)
 mesh_size = 100
 dx = 1.d/double(mesh_size)
 pror1 = dindgen(mesh_size+1)/(mesh_size)
 pror2 = dindgen(mesh_size+2)/(mesh_size) - 1.d/(2.*mesh_size)

 ddbr =  make_array(n_elements(vect[*,0,0])-1,n_elements(vect[0,*,0]),n_elements(vect[0,0,*]),/double)
  for p=0, n_elements(vect[0,0,*])-1 do begin
  for j=0, n_elements(vect[0,*,0])-1 do begin
    dbr = interpol(reform(vect[*,j,p]), pror2, pror1 )
    ddbr[*,j,p] = dbr
  endfor
  endfor
;print,'br interpolated'
  return,ddbr
end
