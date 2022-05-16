function read_spectrum,path
;#; reading the for/eq_csha.blc.for file
 idl_file=path+'for/1/spc3.in'
 dum=file_info(idl_file)

 dissipation=dindgen(2)*0.d
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
;         beta0 = strpos(strtrim(txt(i)),'0')
;         if (total(beta0) gt 0 ) then begin
;          continue
;         endif
;         theta0 = strpos(strtrim(txt(i)),'theta0')
;         if (total(theta0) gt 0 ) then begin
;          continue
;         endif
;         if (total(rr) gt 0 and strcmp(txt(i),'eta0',4) gt 0) then begin
         rr = strpos(strtrim(txt(i)),'mm')
         if (total(rr) gt 0 ) then begin
         dumdum = strmid(txt(i),rr+4)
         lk=strsplit(dumdum,',',/extract)
         mm = fix(lk)
;         print,'long',mm
         endif
        endfor
        for i=0,n_elements(txt)-1 do begin
         red = strpos(strtrim(txt(i)),'nanf')
         if (total(red) gt 0) then begin
;         print,'dove è mue',txt(i)
          dumdum = strmid(txt(i),red+6)
          lk=strsplit(dumdum,',',/extract)
          nanf = fix(lk)
;          print,'bbbb',nanf
         endif
        endfor
        for i=0,n_elements(txt)-1 do begin
         red = strpos(strtrim(txt(i)),'nz')
         if (total(red) gt 0) then begin
;         print,'dove è mue',txt(i)
          dumdum = strmid(txt(i),red+7)
          lk=strsplit(dumdum,',',/extract)
          nz = long(lk)
         endif
        endfor
        for i=0,n_elements(txt)-1 do begin
         red = strpos(strtrim(txt(i)),'nzcon')
         if (total(red) gt 0) then begin
;!         print,'dove è mue',txt(i)
          dumdum = strmid(txt(i),red+7)
          lk=strsplit(dumdum,',',/extract)
          nzcon = long(lk)
         endif
        endfor
      
        aaa = create_struct('mm',mm,'nanf',nanf,'nz',nz,'nzcon',nzcon)
        return,aaa
 endif else begin
;  print,'non si può leggere lo spettro da simulazioni vecchie!'
  restore,path+'dat/spectrum.sav'
  nzcon = nz
  nz(0) = nzcon(0) - 1
  aaa = create_struct('mm',mm,'nanf',nanf,'nz',nz,'nzcon',nzcon)
  return,aaa
;  for_file=path+'for/1/cyl1.inc.for'
;  dum=file_info(for_file)
;  if (dum.exists eq 0) then begin
;   print,'errore, manca il file di configurazione'
;   stop
;  endif else begin
;;   print,for_file
;   nlines = file_lines(for_file)
;   openr,lun,for_file,/get_lun
;   txt=make_array(nlines,/string)
;   readf,lun,txt
;   free_lun,lun
;   strstr='*equilibrium*'
;;   print,txt(where(strmatch(txt,strstr)))
;;   print,where(strmatch(txt,strstr))
;   strstr='*rr = *'
;;   print,txt(where(strmatch(txt,strstr)))
;   pos=where(strmatch(txt,strstr))
;   arr = strpos(txt(pos),'rr')
;;   print,'att',arr
;   dumdum = strmid(txt(pos),arr+5,1)
;;   help,dumdum
;;   print,dumdum
;   RR = float(dumdum)
;   RR = RR * 1.d
;;   print,RR
;   return,dissipation
;  endelse
 endelse
end
