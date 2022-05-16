pro idl_saveb_bin, path, itp

 nsectionsdat = find_nsectionsdat(path)
 xxx = read_settings(path,nsectionsdat)

 print,itp
 help,xxx,/str
 idl_file=path+'dat/ibprofiles.sav'
 print,'restoring...',idl_file
 restore,idl_file
; restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
 
 if (itp lt xxx.min) then begin
  print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
  stop
 endif
 
 ditp=fix(itp-xxx.min) 
 bzt=make_array(n_elements(pdbz[0,*,0]),n_elements(pdbz[0,0,*]),/dcomplex,value=0.)
 btt=make_array(n_elements(pdbt[0,*,0]),n_elements(pdbt[0,0,*]),/dcomplex,value=0.)
 brt=make_array(n_elements(pdbr[0,*,0]),n_elements(pdbr[0,0,*]),/dcomplex,value=0.)
 
 bzt[*,*]=reform(pdbz[ditp,*,*]) 
 btt[*,*]=reform(pdbt[ditp,*,*]) 
 brt[*,*]=reform(pdbr[ditp,*,*]) 
 
;#; first: interpolation of br on the x2 mesh
  d3 = 1
  ddbr =  btt * dcomplex(0.d,0.d)
  help,btt
   for q=0, n_elements(btt[*,0])-1 do begin
   for j=0, n_elements(btt[0,*])-1 do begin
         dbr = interpol(reform(brt[*,j]), pror1, pror2 )
         mn = mnum(j,nz,nanf,d3)
         if (mn[0] eq 0) then begin
          dbr[0] = -dbr[1]
         endif else begin
          if (mn[0] eq 1) then begin
           dbr[0] = dbr[1]
          endif else begin
           dbr[0] = -dbr[1]
          endelse
         endelse
     ddbr[*,j] = dbr
   endfor
   endfor
   help,ddbr

 aaa=size(brt)
 print,'aaa',aaa
 nr = aaa(1)
 modes = aaa(2)


;#; saving 
;#; check directory existence
 dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
 call_procedure,'check_dir',dir 
  
 bin_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/bcyl.bin'
print,'incomputeidllll',bin_file
 print,'writing the magnetic field for itp=',ditp+xxx.min
 openw,lun,bin_file,/get_lun,/f77_unformatted 
 writeu,lun, nr
 writeu, lun, modes
 writeu, lun, ddbr
 writeu, lun, btt
 writeu, lun, bzt
 
 free_lun,lun

end
