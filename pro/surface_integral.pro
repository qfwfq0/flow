function surface_integral,path,vect1,vect2,kind,d3
;#; this function computes the volume integral of a quadratic form like int_dS(A B) 
;input: vect1[*,*],vect2[*,*]: vect1[radius,modes] and vect2[radius,modes]
;vect1, vect2 are complex vectors on the mesh x2

 rr = read_major_radius(path)
 spectrum = read_spectrum(path)
 nz = spectrum.nz
 nanf = spectrum.nanf
 mm = spectrum.mm
 mesh_size = read_mesh(path)
 dx = 1.d/double(mesh_size)
 pror1 = dindgen(mesh_size+1)/(mesh_size)
 pror2 = dindgen(mesh_size+2)/(mesh_size) - 1.d/(2.*mesh_size)
 if (n_elements(vect1[*,0]) ne mesh_size+2) then begin
  print,'errore,mesh radiale di vect1 non è x2',n_elements(vect1[*,0])
  exit
 endif
 if (n_elements(vect2[*,0]) ne mesh_size+2) then begin
  print,'errore,mesh radiale di vect2 non è x2',n_elements(vect2[*,0])
  exit
 endif

fact1 = 2.d * (!Dpi)^2.d * rr(0)
array = make_array(n_elements(vect1[*,0]),/double,value=0.d)

;#; kind = 1 : int_dS(A B) with dS=rd(theta)dz e^r (guscio cilindrico)
case kind of
 1: begin
  for q = 1, n_elements(vect1[*,0]) - 1 do begin ;#; index on radial position
   for p = 0, n_elements(vect1[0,*]) - 1 do begin ;#; Marco.
    mn1 = mnum(p,nz,nanf,d3)
    func = pror1(q-1) * 0.5d * ( (real_part(vect1(q,p))+real_part(vect1(q-1,p))) * ( real_part(vect2(q,p)) + real_part(vect2(q-1,p))) + ( imaginary(vect1(q,p)) + imaginary(vect1(q-1,p))) * ( imaginary(vect2(q,p)) + imaginary(vect2(q-1,p))) )
     
    array(q) = array(q-1) + func
   endfor
  endfor
  array = fact1 * array
 end
 2: begin
  print,'addio'
 end
endcase
;     mn = mnum(j,nz,nanf,d3)
;     cnvl = convo(path,mn[0],mn[1],d3)
;     help,cnvl
;    for p = 0, n_elements(cnvl[0,*])-1 do begin
;     j1 = cnvl[0,p]
;     j2 = cnvl[1,p]
;     mn1 = mnum(j1,nz,nanf,d3)
;     mn2 = mnum(j2,nz,nanf,d3)
;
;     print,'nonlinearconvo0',p,j,j1,j2
;;#; definition of the multiplication factor
;     fact = mlpf(mn1,mn2)
;
;     array = array + fact*vect1[*,j1]*vect2[*,j2] 
;;#; Marco, end of the convolution cycle
;    endfor
    return,array
end
