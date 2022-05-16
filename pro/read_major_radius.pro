function read_major_radius,path
;#; reading the for/eq_csha.blc.for file
 idl_file=path+'for/1/spc3.in'
 dum=file_info(idl_file)

print,idl_file
print,'AAAAA',dum.exists
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
         rr = strpos(txt(i),'rr')
         if (total(rr) gt 0) then begin
          dumdum = strmid(txt(i),rr+4,3)
          RR = float(dumdum)
          break
         endif
        endfor
        RR= RR*1.d
;    print,RR
        if (RR(0) lt 1e-1) then begin
         RR(0)=4.d0
        endif     
        return, RR
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
   strstr='*equilibrium*'
;   print,txt(where(strmatch(txt,strstr)))
;   print,where(strmatch(txt,strstr))
   strstr='*rr = *'
;   print,txt(where(strmatch(txt,strstr)))
   pos=where(strmatch(txt,strstr))
   arr = strpos(txt(pos),'rr')
;   print,'att',arr
   dumdum = strmid(txt(pos),arr+5,1)
;   help,dumdum
;   print,dumdum
   RR = float(dumdum)
   RR = RR * 1.d
   if (RR(0) lt 1e-1) then begin
    RR(0)=4.d0
   endif     
;   print,RR
   return,RR
  endelse
 endelse
end
