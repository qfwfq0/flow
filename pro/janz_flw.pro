;#; this function converts a diad (m,n) to janz
function janz_flw,m,n,my,mm,nanf,nzcon
;#; Marco, modificato il 03/08/2018
;    print,m,n
    j = 0
    if (m eq 0 and n eq 0) then return,j
    if (m eq 0) then begin
     j = nzcon[0] + n
    endif else begin
     for mi=0,m-1 do begin
         j = j + nzcon[mi]
     endfor
     j =  j + n - nanf[m]
    endelse
    return,j
end
