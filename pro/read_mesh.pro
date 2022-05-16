function read_mesh,path
;#; reading the for/eq_csha.blc.for file
 idl_file=path+'for/1/spc2.in'
 dum=file_info(idl_file)

 mesh = 0l
 if (dum.exists eq 1) then begin
;#; first of all I count all the lines in the file
;        print,idl_file
        nlines = file_lines(idl_file)
        openr,lun,idl_file,/get_lun
        txt=make_array(nlines,/string)
        readf,lun,txt
        free_lun,lun

        for i=0,n_elements(txt)-1 do begin
         rr = strpos(strtrim(txt(i)),'lx')
         if (total(rr) gt 0 ) then begin
         dumdum = strmid(txt(i),rr+4)
         lx = long(dumdum)
         mesh = lx
         endif
        endfor
        return, mesh
 endif else begin
   mesh = 100l
  return, mesh
 endelse
end
