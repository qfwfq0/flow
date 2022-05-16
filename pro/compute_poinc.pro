pro compute_poinc, path, d3
;#; I open the idl_quantities.dat file to see when and what I want to compute

;#; first of all I count all the lines in the file
 idl_file=path+'/for/pppoincare.dat'
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
 if (yesorno(0) eq 1) then begin
 endif
  idl_file = path+'dat/ibprofiles.sav'
  cmpt=test_up_to_date2(path,idl_file)
  if (cmpt eq 1) then begin
   mf = 0
   call_procedure,'isavemodeprofiles',path,mf,d3
  endif
  for q = 0, nitp-1 do begin
   ex_file = file_test(path+'dat/itp/'+itp(q)+'/bcyl.bin')   
   if (ex_file ne 1) then begin
    print,'need to create the magnetic field file for itp=',itp(q)
;#; write the magnetic field in unformatted files
    call_procedure,'idl_saveb_bin',path,itp(q)
;    call_procedure,'ppoinc',path,itp(q)
   endif
  endfor

end
