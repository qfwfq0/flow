;#; this function converts janz to a diad (m,n)
function mnum_v2,j,minput,nz,nanf,d3,tok
 if (tok eq 0) then begin ;#; begin RFP case
  if (d3 eq 1) then begin ;#; begin RFP 3D case
;#; compute the cumulated nz array, cnz
     if (j eq 0) then begin
      return,[0,0]
     endif
     cnz = nz
     for q=0, n_elements(nz)-1 do begin
      if (q eq 0) then begin
       cnz(q) = nz(0) ;- 1
      endif else begin
       cnz(q) = cnz(q-1) + nz(q)
      endelse
     endfor
     mm = max(where(j gt cnz)) + 1
     if (mm eq -1) then mm = 0
     if (mm eq 0) then begin
      nn = nanf(mm) + j - 1
     endif else begin
      nn = nanf(mm) + ( j - 1 - cnz( max( where(j gt cnz) ) ) )
     endelse
     return,[mm,nn]
   endif else begin ;#; 2D RFP case
     if (j eq 0) then begin
      return,[0,0]
     endif else begin
      return, [j,nanf(j)]
     endelse
   endelse
 endif else begin ;#; begin tokamak case
  if (d3 eq 1) then begin
   print,'code this mnum_v2.pro'
   stop
  endif else begin ;#; 2d tokamak case
   if (j eq 0) then begin
    return,[0,0]
   endif else begin
    nzne0 = where(nz ne 0)  
;    print,'test',minput
;    print,'test2',nzne0
    return, [minput(nzne0(j)),nanf(nzne0(j))]
   endelse
  endelse
 endelse
end
