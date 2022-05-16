function read_dissipation,path
;#; reading the for/eq_csha.blc.for file

 idl_file=path+'for/1/spc3.in'
 dum1=file_info(idl_file)
;#; Marco 29 agosto, cerco la dissipazione sui file for/3/ perché in certi casi inizio le simulazioni con dissipazioni amggiori per non farle crashare
 dum = make_array(3,/integer,value=-3)
 idl_files = make_array(3,/string)
 idl_files[2]=path+'for/3/spc3.in'
 dum3=file_info(idl_files[2])
 idl_files[1]=path+'for/2/spc3.in'
 dum2=file_info(idl_files[1])
 idl_files[0]=path+'for/1/spc3.in'
 dum1=file_info(idl_files[0])
 dum[0]=dum1.exists
 dum[1]=dum2.exists
 dum[2]=dum3.exists
 qqq=max(where(dum eq 1))
 idl_file = idl_files(qqq)
  


 dissipation=dindgen(2)*0.d
 if (dum1.exists eq 1) then begin
;#; first of all I count all the lines in the file
;        print,idl_file
        nlines = file_lines(idl_file)
        openr,lun,idl_file,/get_lun
        txt=make_array(nlines,/string)
        readf,lun,txt
        free_lun,lun

;#; find RR
        for i=0,n_elements(txt)-1 do begin
         beta0 = strpos(strtrim(txt(i)),'beta0')
         if (total(beta0) gt 0 ) then begin
          continue
         endif
         theta0 = strpos(strtrim(txt(i)),'theta0')
         if (total(theta0) gt 0 ) then begin
          continue
         endif
;         if (total(rr) gt 0 and strcmp(txt(i),'eta0',4) gt 0) then begin
         rr = strpos(strtrim(txt(i)),'eta0')
         if (total(rr) gt 0 ) then begin
         dumdum = strmid(txt(i),rr+6)
         eta0 = float(dumdum)
;         help,eta0
        dissipation(0) = eta0
;          break
         endif
        endfor
        eta0= eta0*1.d
        for i=0,n_elements(txt)-1 do begin
         mm = strpos(strtrim(txt(i)),'mue')
         if (total(mm) gt 0) then begin
;         print,'dove è mue',txt(i)
          dumdum = strmid(txt(i),mm+5)
          mu0 = float(dumdum)
;          help,mu0
;          break
         endif
        endfor
;#; Cerco alet
        for i=0,n_elements(txt)-1 do begin
         mm = strpos(strtrim(txt(i)),'alet')
         if (total(mm) gt 0) then begin
          dumdum = strmid(txt(i),mm+6)
          alet = float(dumdum)
;          print,'alet',alet
;          break
         endif
        endfor
;#; Cerco beet
        for i=0,n_elements(txt)-1 do begin
         mm = strpos(strtrim(txt(i)),'beet')
         if (total(mm) gt 0) then begin
          dumdum = strmid(txt(i),mm+6)
          beet = float(dumdum)
;          print,'beet',beet
;          break
         endif
        endfor
;#; Cerco gaet
        for i=0,n_elements(txt)-1 do begin
         mm = strpos(strtrim(txt(i)),'gaet')
         if (total(mm) gt 0) then begin
          dumdum = strmid(txt(i),mm+6)
          gaet = float(dumdum)
;          print,'gaet',gaet
;          break
         endif
        endfor
;#; Cerco almu
        for i=0,n_elements(txt)-1 do begin
         mm = strpos(strtrim(txt(i)),'almu')
         if (total(mm) gt 0) then begin
          dumdum = strmid(txt(i),mm+6)
          almu = float(dumdum)
;          print,'almu',almu
;          break
         endif
        endfor
;#; Cerco bemu
        for i=0,n_elements(txt)-1 do begin
         mm = strpos(strtrim(txt(i)),'bemu')
         if (total(mm) gt 0) then begin
          dumdum = strmid(txt(i),mm+6)
          bemu = float(dumdum)
;          print,'bemu',bemu
;          break
         endif
        endfor

;#; return quantities
        mu0= mu0*1.d
        alet = alet*1.d
        beet = beet*1.d
        gaet = gaet*1.d
        almu = almu*1.d
        bemu = bemu*1.d
        dissipation(1) = mu0
        aaa = create_struct('eta',eta0,'nu',mu0,'alet',alet,'beet',beet,'gaet',gaet,'almu',almu,'bemu',bemu)
        return, aaa
 endif else begin

  for_file=path+'for/1/cyl1.inc.for'
  dum=file_info(for_file)
  if (dum.exists eq 0) then begin
   print,'errore, manca il file di configurazione'
   stop
  endif else begin
;   print,for_file
   nlines = file_lines(for_file)
   openr,lun,for_file,/get_lun
   txt=make_array(nlines,/string)
   readf,lun,txt
   free_lun,lun
   print,'ciaoooo'
   strstr='*ETA0*'
   qqq=where(strmatch(txt,strstr))
   eta0 = strmid(txt(qqq),17,6)
   eta0=float(eta0)
   strstr='*MUE*'
   qqq=where(strmatch(txt,strstr))
   nu0 = strmid(txt(qqq(-1)),17,6)
   nu0=float(nu0)
   aaa = create_struct('eta',eta0,'nu',nu0)
   return, aaa

;   for i=0,n_elements(txt)-1 do begin
;
;    beta0 = strpos(strtrim(txt(i)),'beta0')
;    if (total(beta0) gt 0 ) then begin
;     continue
;    endif
;    theta0 = strpos(strtrim(txt(i)),'theta0')
;    if (total(theta0) gt 0 ) then begin
;     continue
;    endif
;   strstr='*ETA0*'
;   qqq=where(strmatch(txt,strstr)))
;    rr = strpos(strtrim(txt(i)),'ETA0')
;    if (total(rr) gt 0 ) then begin
;     dumdum = strmid(txt(i),rr+6)
;     eta0 = float(dumdum)
;         help,eta0
;     dissipation(0) = eta0
;;          break
;    endif
;   print,'a.. ',txt(i)
;   strstr='*ETA0*'
;   pos=where(strmatch(txt(i),strstr))
;   if (total(pos) gt 0) then begin
;   arr = strpos(txt(pos),'ETA0')
;   print,'att',txt(pos)
;   dumdum = strmid(txt(pos),arr+5,1)
;   endif else begin
;    continue
;   endelse
;   help,dumdum
;   print,dumdum


  endelse
 endelse
end
