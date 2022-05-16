pro read_spectrum, path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d

 my = 0l
 i2d = 0l
 m2d = 0l
 n2d = 0l
 spectrumfile = path+'dat/spectrum.dat'
 openr, lun, spectrumfile, /get_lun, /f77_unformatted
  readu, lun, my
  mm = l64indgen(my+1)
  nanf = l64indgen(my+1)
  nz = l64indgen(my+1)
  readu, lun, mm 
  readu, lun, nanf 
  readu, lun, nz 
  nzcon = nz
  nzcon(0) = nz(0)+1
  if (d3 eq 0) then begin
   if (eof(lun) eq 1) then begin
    i2d=1
    m2d=1
    n2d=1
   endif else begin
   readu, lun, i2d
   readu, lun, m2d
   readu, lun, n2d
   endelse
  endif
 free_lun, lun
end
