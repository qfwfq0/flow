;#; this function converts a diad (m,n) to janz
function janz2,m,n,my,mm,nanf,nz
;    print,m,n
    j = 0
    if (m eq 0 and n eq 0) then return,j
    if (m eq 0) then begin
     j = nz[0] + n
    endif else begin
     for mi=0,m-1 do begin
         j = j + nz[mi]
     endfor
     j =  j + n - nanf[m]
    endelse
    return,j
end
