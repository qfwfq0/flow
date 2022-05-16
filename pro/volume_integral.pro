	function volume_integral,path,vect1,vect2,kind,d3,tok
	;#; this function computes the volume integral of a quadratic form like int_dV(A B) 
	;input: vect1[*,*],vect2[*,*]: vect1[radius,modes] and vect2[radius,modes]
	;vect1, vect2 are complex vectors on the mesh x2

	;#; kind = 1 : int_dV(A dzB)
	;#; kind = 2 : int_dV(A d(theta)B)
	;#; kind = 3 : int_dV(A B)
	;#; kind = 4 : int_dV(A/r d(theta)B)

	 rr = read_major_radius(path)
	; spectrum = read_spectrum(path)
	; nz = spectrum.nz
	; nzcon = spectrum.nzcon
	; nanf = spectrum.nanf
	; mm = spectrum.mm
	 result = file_test(path+'dat/spectrum.sav')
	 if (result ne 1) then begin
          call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
	 endif else begin
	  restore,path+'dat/spectrum.sav'
          if (nzcon eq !null) then begin
           nzcon = nz
           nzcon(0) = 1
          endif
	 endelse
	 mesh_size = read_mesh(path)
	 dx = 1.d/double(mesh_size)
	 pror1 = dindgen(mesh_size+1)/(mesh_size)
	 pror2 = dindgen(mesh_size+2)/(mesh_size) - 1.d/(2.*mesh_size)
	 if (n_elements(vect1[*,0]) ne mesh_size+2) then begin
	  print,'errore,mesh radiale di vect1 non è x2',n_elements(vect1[*,0])
	  stop
	 endif
	 if (n_elements(vect2[*,0]) ne mesh_size+2) then begin
	  print,'errore,mesh radiale di vect2 non è x2',n_elements(vect2[*,0])
	  stop
	 endif

	;#; fattore comune a tutti gli integrali di volume in dtheta*dz. In più c'è già un fattore 1/2 che serve a tenere conto degli integrali di tipo cos^2 o sin^2 (che sono gli unici che sopravvivono quando integro sul volume).
	;#; REMARK: int_dz(cos^2(z))=pi
	fact1 = 2.d * (!Dpi)^2.d * rr(0)
	;array = make_array(n_elements(vect1[*,0]),/double,value=0.d)
	array = make_array(mesh_size+1,/double,value=0.d)
	func = 0.d

	;#; kind = 1 : int_dV(A dzB)
	case kind of
	 1: begin
	;  for q = 1, n_elements(vect1[*,0]) - 1 do begin ;#; index on radial position
	  for q = 0, mesh_size do begin ;#; index on radial position
	   dum = 0.d
	   func = 0.d
	   for p = 1, n_elements(vect1[0,*]) - 1 do begin ;#; Marco. Il ciclo sui modi parte da p=1 perché la derivata in z del modo 00 vale zero e quindi non contribuisce all'integrale
;	    mn1 = mnum(p,nzcon,nanf,d3)
	    mn1 = mnum_v2(p,mm,nzcon,nanf,d3,tok)
	;#; il fattore 1/4 viene dalla media di vect1*vect2
	    func = - mn1(1) * pror1(q) * 0.25d * ( (real_part(vect1(q+1,p))+real_part(vect1(q,p))) * ( imaginary(vect2(q+1,p)) + imaginary(vect2(q,p))) - ( imaginary(vect1(q+1,p)) + imaginary(vect1(q,p))) * ( real_part(vect2(q+1,p)) + real_part(vect2(q,p))) )
	;    print,mn1(1),( (real_part(vect1(q+1,p))+real_part(vect1(q,p))) * ( imaginary(vect2(q+1,p)) + imaginary(vect2(q,p))) - ( imaginary(vect1(q+1,p)) + imaginary(vect1(q,p))) * ( real_part(vect2(q+1,p)) + real_part(vect2(q,p))) )
	    dum = dum + func * dx * fact1
	   endfor
	    if (q eq 0) then begin
	     array(q) =  dum
	    endif else begin
	     array(q) = array(q-1) + dum
	    endelse
	  endfor
	;  print,'inizio',array
	 end
	;#; kind = 2 : int_dV(A d(theta)B)
	 2: begin
	;  for q = 1, n_elements(vect1[*,0]) - 1 do begin ;#; index on radial position
	  for q = 0, mesh_size do begin ;#; index on radial position
	   dum = 0.d
	   func = 0.d
	   for p = nz(0)+1, n_elements(vect1[0,*]) - 1 do begin ;#; Marco. Il ciclo sui modi parte da p=nz(0) perché la derivata in theta del modo 0n vale zero e quindi non contribuisce all'integrale
;	    mn1 = mnum(p,nzcon,nanf,d3)
	    mn1 = mnum_v2(p,mm,nzcon,nanf,d3,tok)
	;#; il fattore 1/4 viene dalla media di vect1*vect2
	    func =  - mn1(0) * pror1(q-1) * 0.25d * ( (real_part(vect1(q+1,p))+real_part(vect1(q,p))) * ( imaginary(vect2(q+1,p)) + imaginary(vect2(q,p))) - ( imaginary(vect1(q+1,p)) + imaginary(vect1(q,p))) * ( real_part(vect2(q+1,p)) + real_part(vect2(q,p))) )
	    dum = dum + func * fact1 * dx 
	   endfor
	    if (q eq 0) then begin
	     array(q) =  dum
	    endif else begin
	     array(q) = array(q-1) + dum
	    endelse
	  endfor
	 end
	;#; kind = 3 : int_dV(A B)
	 3: begin
	;  for q = 1, n_elements(vect1[*,0]) - 1 do begin ;#; index on radial position
	  for q = 0, mesh_size do begin ;#; index on radial position
	   dum = 0.d
	   for p = 0, n_elements(vect1[0,*]) - 1 do begin ;#; ciclo sui modi
	    mn1 = mnum_v2(p,mm,nzcon,nanf,d3,tok)
	;#; il fattore 1/4 viene dalla media di vect1*vect2
    func = pror1(q-1) * 0.25d * ( (real_part(vect1(q+1,p))+real_part(vect1(q,p))) * ( real_part(vect2(q+1,p)) + real_part(vect2(q,p))) + ( imaginary(vect1(q+1,p)) + imaginary(vect1(q,p))) * ( imaginary(vect2(q+1,p)) + imaginary(vect2(q,p))) )
     
     if (p eq 0) then begin
;#; Marco, aprile 2018: avevo messo fact1*2 nei modi diversi da (0,0), mentre è il contrario. vedi pagina marcata AA nel quaderno Aix17
      dum = dum + func * fact1 * dx * 2.d
     endif else begin
      dum = dum + func * fact1 * dx ;* 2.d
     endelse 
   endfor

   if (q eq 0) then begin
    array(q) =  dum
   endif else begin
    array(q) = array(q-1) + dum
   endelse
;   print,'q',array(q),dum
  endfor
 end
;#; kind = 4 : int_dV(A/r d(theta)B)
 4: begin
;  print,'caso4'
;  for q = 1, n_elements(vect1[*,0]) - 1 do begin ;#; index on radial position
  for q = 0, mesh_size do begin ;#; index on radial position
   dum = 0.d
   func = 0.d
   for p = nz(0)+1, n_elements(vect1[0,*]) - 1 do begin ;#; Marco. Il ciclo sui modi parte da p=nz(0) perché la derivata in theta del modo 0n vale zero e quindi non contribuisce all'integrale
;    mn1 = mnum(p,nzcon,nanf,d3)
    mn1 = mnum_v2(p,mm,nzcon,nanf,d3,tok)
	;#; il fattore 1/4 viene dalla media di vect1*vect2
    func =  - mn1(0) * 0.25d * ( (real_part(vect1(q+1,p))+real_part(vect1(q,p))) * ( imaginary(vect2(q+1,p)) + imaginary(vect2(q,p))) + ( imaginary(vect1(q+1,p)) + imaginary(vect1(q,p))) * ( real_part(vect2(q+1,p)) + real_part(vect2(q,p))) )

    dum = dum + func * fact1 * dx
   endfor

   if (q eq 0) then begin
    array(q) =  dum
   endif else begin
    array(q) = array(q-1) + dum 
   endelse

  endfor
 end
endcase

    return,array
end
