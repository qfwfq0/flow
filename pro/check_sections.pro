pro check_sections,path,mf,d3
; scopro il numero di sezioni dat e sav
    nsectionsdat = find_nsectionsdat(path)
    nsectionssav = find_nsectionssav(path)
    lsectionsdat = intarr(nsectionsdat)
    print,'numero sezioni dat: ', nsectionsdat
    print,'numero sezioni sav: ', nsectionssav
    print,'d3',d3
;#; array containg the number of time steps written on the .dat file
    written_itp = ulonarr(nsectionsdat)
;#; array containing the number of time steps computed by the code
    computed_itp = ulonarr(nsectionsdat)
    xyz=read_settings(path,nsectionsdat)
    computed_itp = xyz.computed
    max_itp = xyz.max

;#; here I check if some .sav sections are missing
    if (nsectionsdat gt nsectionssav) then begin
     missing = nsectionsdat - nsectionssav
     print,'missing',missing
     for q=0,missing-1 do begin
      ;#; typically the last dat is missing
      qq = nsectionsdat
      filein = 'specyl_'+strtrim(string(qq,format='(i0)'))+'_all.dat'
      fileout = 'specyl_'+strtrim(string(qq,format='(i0)'))+'_all.sav'
      print,'ireadscalling:      ',filein,' ----> ',fileout,qq,computed_itp(qq-1)
      call_procedure,'ireadsc',path,fileout,filein,qq,qq,mf,d3
     endfor
;#; now the number of sav files is equal to the number of dat files and we can go on
      nsectionssav = nsectionsdat
    endif

    if (nsectionssav gt 0) then begin
     lsections = intarr(nsectionssav)
    endif

;#; here I check wheter the highest specyl_?_all.sav file contains the required number of time steps written in computed_itp(n_elements(computed_itp)-1)
   
;#; slow part!
;        file = path+string(nsectionssav,format=fileformat)
;        print,'éééééééééééééééééé',file
;        restore,file
;        qqq=n_elements(bz[*,0,0])-1
;;        qqq=computed_itp(nsectionssav - 1)
;;        if (qqq eq computed_itp(nsectionssav - 1)) then begin
;        print,'checaz',qqq,max_itp
;        if (qqq eq max_itp) then begin
;         print,'everything is fine',qqq,computed_itp(nsectionssav-1)
;        endif else begin
;         print,'I need to recompute the highest rank actually present specyl_?_all.sav',fileout,computed_itp(nsectionssav - 1)
;         filein = 'specyl_'+strtrim(string(nsectionssav,format='(i0)'))+'_all.dat'
;         fileout = 'specyl_'+strtrim(string(nsectionssav,format='(i0)'))+'_all.sav'
;         call_procedure,'ireadsc',path,fileout,filein,nsectionssav - 1
;        endelse

;#; here I know how many time steps I want to read from each specyl_?_all.dat file, and I can safely call the ireadsc procedure.
        
       for q = nsectionssav + 1, nsectionsdat do begin
        print,q
         filein = 'specyl_'+strtrim(string(q,format='(i0)'))+'_all.dat'
         fileout = 'specyl_'+strtrim(string(q,format='(i0)'))+'_all.sav'
         print,'I need to recompute',fileout,computed_itp(q - 1)
;:        call_procedure,'ireadsc',path,fileout,filein,computed_itp(q - 1)
       endfor

end


;#; parts of the code that should be in pro/read_settings.pro
;        i = 1
;;#; reading the for/settings.for file
;;#; first of all I count all the lines in the file
;        idl_file=path+'for/settings.blc.for'
;        print,idl_file
;        nlines = file_lines(idl_file)
;        openr,lun,idl_file,/get_lun
;        txt=make_array(nlines,/string)
;        readf,lun,txt
;        free_lun,lun
;;#; minitp
;        k0 = strpos(txt(4),'/')
;        dum = txt(4)
;        dumdum = strmid(dum,k0+1)
;        k0 = strpos(dumdum,'/')
;        minitp = fix(strmid(dumdum,0,k0))
;        print,'min itp required=',minitp
;
;;#; maxitp
;        k0 = strpos(txt(5),'/')
;        dum = txt(5)
;        dumdum = strmid(dum,k0+1)
;        k0 = strpos(dumdum,'/')
;        maxitp = fix(strmid(dumdum,0,k0))
;        print,'max itp required=',maxitp
;
;;#; the first specyl_?_all.dat file ends at...
;        k0 = strpos(txt(6),'/')
;        k1 = strpos(txt(6),',')
;        if (k0 eq -1 or k1 eq -1) then begin
;         print,'there must be an error!'
;        endif
;        dum=txt(6)
;        dumdum=strmid(dum,k0+1,k1-k0-1)
;        print,dumdum
;        computed_itp(0)=fix(dumdum)
;        newdum=strmid(dum,k1)
;        counter = 1
;;#; here I have read how many time steps are computed in the first dat_section.
;;#; then we find how many time steps are computed on the others specyl_?_all.dat
;        while (counter le nsectionsdat - 1) do begin
;;#; here I remove the comma
;         comp = strmid(newdum,strpos(newdum,',')+1,strlen(newdum))
;;#; here I read the number
;         readn = strmid(comp,0,strpos(comp,','))
;         if (fix(readn) eq 0) then begin
;          computed_itp(counter) = maxitp - total(computed_itp)
;         endif else begin
;          computed_itp(counter) = fix(readn) - total(computed_itp)
;         endelse
;         counter = counter + 1
;         newdum = strmid(comp,strlen(readn),strlen(newdum))
;;         if (i ge nsectionsdat) then break
;        endwhile
;;        endif
; 
;;#; I check wheter SpeCyl has computed all the time steps that I asked for, through maxitp
;        if (maxitp gt total(computed_itp)) then begin
;         print,'there must be an error, either on for/settings or on the code output.'
;         stop
;        endif
;        print,'ITP computed by SpeCyl',computed_itp


