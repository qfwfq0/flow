function read_settings,path,nsectionsdat
;#; array containg the number of time steps written on the .dat file
    written_itp = ulonarr(nsectionsdat)
;#; array containing the number of time steps computed by the code
    computed_itp = ulonarr(nsectionsdat)

        i = 1
;#; reading the for/settings.for file
;#; first of all I count all the lines in the file
        idl_file=path+'for/settings.blc.for'
        print,idl_file
        nlines = file_lines(idl_file)
        openr,lun,idl_file,/get_lun
        txt=make_array(nlines,/string)
        readf,lun,txt
        free_lun,lun
;#; minitp
        k0 = strpos(txt(4),'/')
        dum = txt(4)
        dumdum = strmid(dum,k0+1)
        k0 = strpos(dumdum,'/')
        minitp = fix(strmid(dumdum,0,k0))
        print,'min itp required=',minitp

;#; maxitp
        k0 = strpos(txt(5),'/')
        dum = txt(5)
        dumdum = strmid(dum,k0+1)
        k0 = strpos(dumdum,'/')
        maxitp = fix(strmid(dumdum,0,k0))
        print,'max itp required=',maxitp

;#; the first specyl_?_all.dat file ends at...
        k0 = strpos(txt(6),'/')
        k1 = strpos(txt(6),',')
        if (k0 eq -1 or k1 eq -1) then begin
         print,'there must be an error!'
        endif
        dum=txt(6)
        dumdum=strmid(dum,k0+1,k1-k0-1)
        print,dumdum
        computed_itp(0)=fix(dumdum)
        newdum=strmid(dum,k1)
        counter = 1
;#; here I have read how many time steps are computed in the first dat_section.
;#; then we find how many time steps are computed on the others specyl_?_all.dat
        while (counter le nsectionsdat - 1) do begin
;#; here I remove the comma
         comp = strmid(newdum,strpos(newdum,',')+1,strlen(newdum))
;#; here I read the number
         readn = strmid(comp,0,strpos(comp,','))
         if (fix(readn) eq 0) then begin
          computed_itp(counter) = maxitp - total(computed_itp)
         endif else begin
          computed_itp(counter) = fix(readn) - total(computed_itp)
         endelse
         counter = counter + 1
         newdum = strmid(comp,strlen(readn),strlen(newdum))
;         if (i ge nsectionsdat) then break
        endwhile
;        endif
 
;#; I check wheter SpeCyl has computed all the time steps that I asked for, through maxitp
        if (maxitp gt total(computed_itp)) then begin
         print,'there must be an error, either on for/settings or on the code output.'
         stop
        endif
        print,'ITP computed by SpeCyl',computed_itp

        cumulated = computed_itp
        for i=1,n_elements(computed_itp)-1 do begin
         cumulated(i) = cumulated(i-1) + computed_itp(i)
        endfor

        aaa = create_struct('min',long(minitp),'max',long(maxitp),'computed',long(computed_itp),'cumulated',long(cumulated))
        return,aaa

end
