pro compute_idl, path, d3
;#; I open the idl_quantities.dat file to see when and what I want to compute

;#; first of all I count all the lines in the file
 idl_file=path+'/for/idl_quantities.dat'
 nlines = file_lines(idl_file)
 openr,lun,idl_file,/get_lun
 txt=make_array(nlines,/string)
 readf,lun,txt
 free_lun,lun
;#; then I remove the commented lines
 cnt = 0
 useful = make_array(nlines,/long,value=-2)
 for q = 0, nlines-1 do begin
  k = strpos(txt(q),'#')
;  print,q,k
  if (k eq -1) then begin
   useful(cnt) = q
   cnt = cnt+1
  endif
 endfor
 text = txt(useful(where(useful ge 0)))
 for i=0,n_elements(text)-1 do begin
;  print,text(i)
 endfor
 
 division = where(text eq '$')
;#; division[0] - 1 = number of itp to be computed 
;#; float(text(division(1)-1)) = chosen helicity
;#; float(division(2) - division(1) -2 = number of available quantities

;#; creation of an array with the chosen times
 nitp = fix(division(0)) - 1
 itp = make_array(nitp,/string)
 for q = 0, nitp - 1 do begin
  itp(q) = text(1+q)
 endfor
;#; creation of the ieli (helicity) quantity
 ieli = text(division(1)-1)
;#; number of available quantities
 nquant = fix(division(2)) - fix(division(1)) - 2
 print,nquant
;#; array containing the available quantites
;#; array containing wheter I want or not the quantity
 quant = make_array(nquant,/string)
 yesorno = make_array(nquant,/integer)

 for q = 0, nquant-1 do begin
  quant[q] = text(division(1) + 2 + q)
  yesorno[q] = fix(text(division(2) + 2 + q))
 endfor

 bcmp=where(yesorno gt 0)
 nbcmp=n_elements(where(yesorno gt 0))
 for q=0, nbcmp-1 do begin
  print,'I need to compute'+quant(bcmp(q))
 endfor

;#; create the directory structure
  for q = 0, nitp-1 do begin
   spawn, 'mkdir -p '+path+'dat/itp/'+itp(q)
  endfor
;#; if yesorno(0)=1 I need to compute [bz,bt,br] 
;#; if yesorno(1)=1 I need to compute [vz,vt,vr] 
;#; if yesorno(2)=1 I need to compute [az,at] 
;#; if yesorno(3)=1 I need to compute [jz,jt,jr] 
;#; if yesorno(4)=1 I need to compute [chi], the helical flux function
;#; if yesorno(5)=1 I need to compute [taujb], the electromagnetic torque
;#; if yesorno(6)=1 I need to compute [flowxy], vector

;#; if yesorno(0)=1 I need to compute [bz,bt,br] 
 if (yesorno(0) eq 1) then begin
  idl_file = path+'dat/ibprofiles.sav'
  call_procedure,'test_up_to_date',path,idl_file,'isavemodeprofiles'
  for q = 0, nitp-1 do begin
   ex_file = file_test(path+'dat/itp/'+itp(q)+'/bcyl.sav')   
   if (ex_file ne 1) then begin
    print,'need to create the magnetic field file for itp=',itp(q)
    call_procedure,'idl_saveb',path,itp(q)
   endif
;#; remember to trate the case exist??ibprofiles? maybe in the makefile
;#; Marco, ask for plot.
    read,in,prompt='Do you want to plot it? Enter 1 if yes,0 if no  '
    if ((in ne 0) and (in ne 1)) then begin
     print,'type 0 (if it is not OK )or 1 (if it is OK)!',in
    endif else begin
     if (in eq 0) then begin  
      exit
     endif else begin
      read,m,prompt='Write the m of the mode'
      print,'m=',m
      read,n,prompt='Write the n of the mode'
      print,'n=',n
      call_procedure,'idl_pltb',path,itp(q),m,n,d3
      
     endelse
    endelse
  endfor
 endif

;#; if yesorno(1)=1 I need to compute [vz,vt,vr] 
 if (yesorno(1) eq 1) then begin
  idl_file = path+'dat/ivprofiles.sav'
  call_procedure,'test_up_to_date',path,idl_file,'isavemodeprofiles'
  for q = 0, nitp-1 do begin
   ex_file = file_test(path+'dat/itp/'+itp(q)+'/vcyl.sav')   
   if (ex_file ne 1) then begin
    print,'need to create the velocity field file for itp=',itp(q)
    call_procedure,'idl_savev',path,itp(q)
   endif
  endfor
 endif

;#; if yesorno(2)=1 I need to compute [az,at] 
 if (yesorno(2) eq 1) then begin
  for q = 0, nitp-1 do begin
   ex_file = file_test(path+'dat/itp/'+itp(q)+'/acyl.sav')   
   if (ex_file ne 1) then begin
    print,'need to create the vector potential field file for itp=',itp(q)
    call_procedure,'idl_savea',path,itp(q),d3
   endif
  endfor
 endif

;#; if yesorno(3)=1 I need to compute [jz,jt,jr] 
 if (yesorno(3) eq 1) then begin
  for q = 0, nitp-1 do begin
   ex_file = file_test(path+'dat/itp/'+itp(q)+'/jcyl.sav')   
   if (ex_file ne 1) then begin
    print,'need to create the current density field file for itp=',itp(q)
    call_procedure,'idl_savej',path,itp(q)
   endif
  endfor
 endif

;#; if yesorno(4)=1 I need to compute [chi], the helical flux function
 if (yesorno(4) ge 1) then begin
;#; check wheter the vector potential, at the selected time,  is computed
  print,itp
  for q = 0, nitp-1 do begin
   ex_file = file_test(path+'dat/itp/'+itp(q)+'/acyl.sav')   
  print,ex_file
   if (ex_file ne 1) then begin
    print,'need to create the vector potential field file for itp= ',itp(q)
    call_procedure,'idl_savea',path,itp(q),d3
   endif
   ex_file = file_test(path+'dat/itp/'+itp(q)+'/chi'+strtrim(string(ieli,format='(i0)'))+'.sav')
   if (ex_file ne 1) then begin
    print,'need to create the helical flux function for itp= ',itp(q),'  3D=',d3
    nth = 128
    nzz = 512
;#; if d3=1 then I need chi on a 3D grid, otherwise on a 2D slice (either poloidal or toroidal)
;    d3 = 1
    call_procedure,'idl_savechi',path,itp(q),ieli,d3,nth,nzz
   endif
  endfor
  
 endif

;#; if yesorno(5)=1 I need to compute [taujb], the electromagnetic torque
 if (yesorno(5) eq 1) then begin
  for q = 0, nitp-1 do begin
   print, 'tau',itp(q)
   call_procedure,'tau',path,d3,itp(q)
  endfor
 endif

;#; if yesorno(6)=1 I need to compute [flowxy], vector of m=1 n=hel flow
 if (yesorno(6) eq 1) then begin
  read,n,prompt='Write the helicity of the m=1 mode '
  print,'n=',n
  for q = 0, nitp-1 do begin
   print, 'flowxy',itp(q)
   call_procedure,'flowxy',path,d3,itp(q),n
   call_procedure,'plot_flw',path,d3,itp(q),n
  endfor
 endif
end
