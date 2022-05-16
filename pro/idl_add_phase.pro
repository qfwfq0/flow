function idl_add_phase,br,modif,delta,min,max

  dimension = size(br)
  print,dimension
  if (dimension[0] eq 2) then begin ;#; 2D sim 
;#; classify the angle
   case 1 of
    (delta gt 0.) and (delta lt !pi): sig_a = -1.
    (delta ge !pi) and (delta lt 2.*!pi): sig_a = 1.
   endcase
   case 1 of
    (delta gt 0.) and (delta lt !pi/2.): sig_b = +1.
    (delta ge !pi/2.) and (delta lt 3.*!pi/2.): sig_b = -1.
    (delta ge 3.*!pi/2.) and (delta lt 2.*!pi): sig_b = +1.
   endcase

   add_ph = strcmp(modif,'r')
   if (add_ph eq 1) then begin
;#; here I modify the radial component of the field (coming from the old version of SpeCyl)
;#; Marco, add a phase to the field: I add a (radially) constant phase the component of the field. 
;#; I'll re-write the fields as:
 ;#; br = alfa_r*im(brold) + i*(1-alfa_r)*im(brold)
 ;#; bt,z = alfa_tz*re(bt,zold) + i*(1-alfa_tz)*re(bt,zold)

   for q=0, n_elements(br[*,0])-1 do begin
    for p = min, max do begin
    if (imaginary(br[q,p]) ge 0.) then begin
     alfa = sig_a * sqrt( (imaginary(br[q,p]) * tan(delta))^2. / (1.+tan(delta)^2.) )
     beta = sig_b * (-1) * alfa /tan(delta)
     br[q,p] = dcomplex(alfa,beta)
    endif else begin
     alfa = -sig_a * sqrt( (imaginary(br[q,p]) * tan(delta))^2. / (1.+tan(delta)^2.) )
     beta = -sig_b * (-1) * alfa /tan(delta)
     br[q,p] = dcomplex(alfa,beta)
    endelse
    endfor
   endfor

   return, br

   endif else begin
   bdd_ph = strcmp(modif,'t')
    if (bdd_ph eq 1) then begin
 ;#; here I modify the poloidal component of the field (coming from the old version of SpeCyl)
 ;#; Marco, add a phase to the field: I add a (radially) constant phase the component of the field. 
 ;#; I'll re-write the fields as:
  ;#; br = alfa_r*im(brold) + i*(1-alfa_r)*im(brold)
  ;#; bt,z = alfa_tz*re(bt,zold) + i*(1-alfa_tz)*re(bt,zold)
   for q=0, n_elements(br[*,0])-1 do begin
    for p = min, max do begin
    if (real_part(br[q,p]) ge 0.) then begin
     alfa = sig_b * sqrt( (real_part(br[q,p]))^2. / (1.+tan(delta)^2.) )
     beta = -sig_a * alfa /tan(delta)
     br[q,p] = dcomplex(alfa,beta)
    endif else begin
     alfa = -sig_b * sqrt( (real_part(br[q,p]))^2. / (1.+tan(delta)^2.) )
     beta = sig_a * alfa /tan(delta)
     br[q,p] = dcomplex(alfa,beta)
    endelse
    endfor
   endfor

   return, br

    endif else begin
 ;#; here I modify the axial component of the field (coming from the old version of SpeCyl)
 ;#; Marco, add a phase to the field: I add a (radially) constant phase the component of the field. 
 ;#; I'll re-write the fields as:
  ;#; br = alfa_r*im(brold) + i*(1-alfa_r)*im(brold)
  ;#; bt,z = alfa_tz*re(bt,zold) + i*(1-alfa_tz)*re(bt,zold)
   for q=0, n_elements(br[*,0])-1 do begin
    for p = min, max do begin
    if (real_part(br[q,p]) ge 0.) then begin
     alfa = sig_b * sqrt( (real_part(br[q,p]))^2. / (1.+tan(delta)^2.) )
     beta = -sig_a * alfa /tan(delta)
     br[q,p] = dcomplex(alfa,beta)
    endif else begin
     alfa = -sig_b * sqrt( (real_part(br[q,p]))^2. / (1.+tan(delta)^2.) )
     beta = sig_a * alfa /tan(delta)
     br[q,p] = dcomplex(alfa,beta)
    endelse
    endfor
   endfor

   return, br


    endelse
   endelse
  endif else begin ;#; 3D sim
  endelse

end
