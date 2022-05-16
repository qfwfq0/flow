function read_momsour,path
;#; reading the for/eq_csha.blc.for file
 idl_file=path+'for/1/spc3.in'
 dum=file_info(idl_file)

 momsour=dindgen(3)*0.d
 if (dum.exists eq 1) then begin
;#; first of all I count all the lines in the file
;        print,idl_file
        nlines = file_lines(idl_file)
        openr,lun,idl_file,/get_lun
        txt=make_array(nlines,/string)
        readf,lun,txt
        free_lun,lun

;#; find RR
        for i=0,n_elements(txt)-1 do begin
         m0_momsour = strpos(strtrim(txt(i)),'m0_momsour')
         if (total(m0_momsour) gt 0 ) then begin
         beg = m0_momsour + 10 + 2
         dumdum = strmid(txt(i),beg)
         re=strsplit(dumdum,/extract)
         aaa = double(re)
         endif
        endfor

endif else begin
      print,'non esiste il file spc3.in con scritto dentro il valore della sorgente di momento'
      stop
endelse

   return, aaa

end
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


