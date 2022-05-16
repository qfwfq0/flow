function janz,m,n,my,mm,nanf,nz
    j = 0
    if (m eq 0 and n eq 0) then return,j
    for mi=0,m-1 do begin
        j = j + nz[mi]
    endfor
    j =  j + 1 + n - nanf[m]
    return,j
end


pro isaveeq,path

    path0 = path
    path = path
    RR = read_major_radius(path)

;    restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
    jmin = 0
    jmax = 0
    print,'estremi j che salvo: ',jmin,jmax
    print,'(ricorda che gli estremi j vanno cambiati a mano!)'

;#; controllo che i specyl_?_all.sav siano a posto.

; scopro il numero di sezioni dat
    nsectionsdat = find_nsectionsdat(path)
    lsectionsdat = intarr(nsectionsdat)
    print,'numero sezioni dat: ', nsectionsdat
    call_procedure,'check_sections',path

;;; scopro il numero di sezioni sav
;    nsectionssav = 1
;    fileformat = '("dat/specyl_",I0,"_all.sav")'
;    while (file_test(path+string(nsectionssav,format=fileformat))) do begin
;        print,'guardo se esiste la sezione ',nsectionssav
;;        print,'file: ',path+string(nsections,format=fileformat)
;;        print,'risultato: ',file_test(path+string(nsections,format=fileformat))
;        nsectionssav = nsectionssav + 1
;    endwhile
;;;    print,'nsections',nsections
;    nsectionssav = nsectionssav - 1
;    print,'numero sezioni sav: ', nsectionssav

    nt = ulonarr(nsectionssav)
    for i = 1,nsectionssav do begin
        file = path+string(i,format=fileformat)
        restore,file
    
        dr = 1.d/lx
        pror1 = dindgen(lx+1)*dr
        pror2 = (-0.5 + dindgen(lx+2))*dr
        nt[i-1] = n_elements(t)
        print,file,nt
        help,bz
        help,bt
        help,br

        pbz = reform(2.*(bz[0:nt[i-1]-1,*,jmin:jmax]))
        pbt = reform(2.*(bt[0:nt[i-1]-1,*,jmin:jmax]))
        pbr = reform(2.*(br[0:nt[i-1]-1,*,jmin:jmax]))
        pvz = reform(2.*(vz[0:nt[i-1]-1,*,jmin:jmax]))
        pvt = reform(2.*(vt[0:nt[i-1]-1,*,jmin:jmax]))
        pvr = reform(2.*(vr[0:nt[i-1]-1,*,jmin:jmax]))
        pt = reform(t(0:nt[i-1]-1))
        help,pbz
        help,pvr

;#; I remove the 2. factor from the equilibrium harmonic
        if (jmin eq 0) then begin
            pbz[*,*,0] = 0.5*pbz[*,*,0]
            pbt[*,*,0] = 0.5*pbt[*,*,0]
            pbr[*,*,0] = 0.5*pbr[*,*,0]
            pvz[*,*,0] = 0.5*pvz[*,*,0]
            pvt[*,*,0] = 0.5*pvt[*,*,0]
            pvr[*,*,0] = 0.5*pvr[*,*,0]
        endif

        if (i eq 1) then begin
            pdbz = temporary(pbz)
            pdbt = temporary(pbt)
            pdbr = temporary(pbr)
            pdvz = temporary(pvz)
            pdvt = temporary(pvt)
            pdvr = temporary(pvr)
            pdt = temporary(pt)

            print,'memory1',memory(/current)
        endif else begin

            pdbz=[pdbz,temporary(pbz[1:*,*,*])]
            pdbt=[pdbt,temporary(pbt[1:*,*,*])]
            pdbr=[pdbr,temporary(pbr[1:*,*,*])]
            pdvz=[pdvz,temporary(pvz[1:*,*,*])]
            pdvt=[pdvt,temporary(pvt[1:*,*,*])]
            pdvr=[pdvr,temporary(pvr[1:*,*,*])]
            pdt = [pdt,temporary(pt[1:*])]
            print,'memory2',memory(/current)

        endelse

        bz = 0
        bt = 0
        br = 0
        vz = 0
        vt = 0
        vr = 0
        pt = 0
;#; end of the main for cycle
    endfor
        

   print,'Sto salvando il file.. attendere.. ib00profiles'

   pdbz=real_part(pdbz)
   pdbt=real_part(pdbt)
   save_file = path+'dat/ib00profiles.sav'
   save,pdt,pror2,pdbz,pdbt,filename=save_file

   pdvz=pdbz
   pdvt=pdbt
   pdvr=pdbr
   save_file = path+'dat/iv00profiles.sav'
   save,pdt,pror2,pdvr,pdvt,pdvz,filename=save_file

end
