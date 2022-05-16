pro write_to_fortran,path,rad1,rad2,vecr,vect,vecz,itp

 nr2 = n_elements(vecr[*,0])
 nr1 = n_elements(vecz[*,0])
 nmodes = n_elements(vecz[0,*])
 help,itp

 openw, lun, path+'dat/itp/'+strtrim(string(itp,'(i0)'))+'/jprof_'+strtrim(string(itp,'(i0)'))+'.dat', /get_lun, /f77_unformatted
 itp = long64(itp)
 writeu,lun,itp
 nmodes = long64(nmodes)
 writeu,lun,nmodes
 print,'nmodes=',nmodes
 nr1 = long64(nr1)
 writeu,lun,nr1
 nr2 = long64(nr2)
 writeu,lun,nr2

 for k=0,nmodes-1 do begin
  for p=0,nr2-1 do begin
   dum = vecr[p,k]
   writeu,lun,dum
  endfor
 endfor
 for k=0,nmodes-1 do begin
  for p=0,nr1-1 do begin
   cum = vect[p,k]
   eum = vecz[p,k]
   writeu,lun,cum
   writeu,lun,eum
  endfor
 endfor

 free_lun,lun
 save,filename=path+'dat/itp/'+strtrim(string(itp,'(i0)'))+'/jprof_'+strtrim(string(itp,'(i0)'))+'.sav',vecr,vect,vecz
end
