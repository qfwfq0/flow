;#; this function converts janz to a diad (m,n)
function mnum,j,nz,nanf,d3
 if (d3 eq 1) then begin
;#; compute the cumulated nz array, cnz
    if (j eq 0) then begin
     return,[0,0]
    endif
    cnz = nz
    for q=0, n_elements(nz)-1 do begin
     if (q eq 0) then begin
      cnz(q) = nz(0) - 1
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
  endif else begin
;#; 2D case
    if (j eq 0) then begin
     return,[0,0]
    endif else begin
     return, [j,nanf(j)]
    endelse
  endelse
end
