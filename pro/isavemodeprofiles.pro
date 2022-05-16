;function janz,m,n,my,mm,nanf,nz
;    j = 0
;    if (m eq 0 and n eq 0) then return,j
;    for mi=0,m-1 do begin
;        j = j + nz[mi]
;    endfor
;    j =  j + 1 + n - nanf[m]
;    return,j
;end


pro isavemodeprofiles,path,mf,d3

    path0 = path
    path = path
    RR = read_major_radius(path)

;    restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
    jmin = 0
;    jmax = 1
;    jmax = 4
;    jmax = 90
;    jmax = 225
   if (d3 eq 1) then begin
;    jmax = 90
    jmax = 135
   endif else begin
    jmax = 10
   endelse
    print,'estremi j che salvo: ',jmin,jmax
    print,'(ricorda che gli estremi j vanno cambiati a mano!)'

;#; controllo che i specyl_?_all.sav siano a posto.

; scopro il numero di sezioni dat e sav
    nsectionsdat = find_nsectionsdat(path)
    nsectionssav = find_nsectionssav(path)
    lsectionsdat = intarr(nsectionsdat)
    print,'numero sezioni dat: ', nsectionsdat
    print,'numero sezioni sav: ', nsectionssav

    xxx=read_settings(path,nsectionsdat)
    minitp = xxx.min
    maxitp = xxx.max
    computed_itp = xxx.computed
    cumulated = xxx.cumulated

;#; Marco, 03/05 commento questa parte vecchia
    print,'controlla che tutte le sezioni siano in ordine'
;    call_procedure,'check_sections',path,mf,d3

;#; given the minimum itp this is the starting section
    min_sec = min(where(minitp lt cumulated) + 1)
;#; given the maximum itp this is the ending section
    max_sec = min(where(maxitp le cumulated) + 1)
;#; Marco, marzo 2018. ridefinisco max_sec
    max_sec = nsectionssav

    print,'min_sec',min_sec
    print,'max_Sec',max_sec
;#; main cycle on the sections
; lettura dal file all.sav
;#; devo leggere da min_sec a partire da minitp-cumulated(min_sec-1)
    if (min_sec ge 2) then begin
     st_time = minitp - cumulated(min_sec - 2)
     end_time = maxitp - cumulated(max_sec - 2)
    endif else begin
     st_time = minitp
     end_time = maxitp
    endelse
    print,'starting_time',st_time,end_time

;    read,in,prompt='Are you sure you want to recompute i*modeprofiles.sav? Enter 1 if yes,0 if no  '
    in = 1
    if ((in ne 0) and (in ne 1)) then begin
     print,'type 0 (if it is not OK )or 1 (if it is OK)!',in
    endif else begin
     if ( in eq 1 or in eq "y") then begin
      print,'yes!'
     endif else begin
      print,'no!'
      stop
     endelse
    endelse

    nt = ulonarr(nsectionsdat)
    fileformat = '("dat/specyl_",I0,"_all.sav")'
;    for i = 1,nsectionssav do begin

    print,'inizia il ciclo',min_sec,max_sec
    print, ' '
    for i = min_sec[0],max_sec[0] do begin
        file = path+string(i,format=fileformat)
        restore,file
        jmax = n_elements(bz[0,0,*]) - 1
    
        dr = 1.d/lx
        pror1 = dindgen(lx+1)*dr
        pror2 = (-0.5 + dindgen(lx+2))*dr
        nt[i-1] = n_elements(t)

        if (i eq min_sec) then begin
         mint = st_time
         maxt = n_elements(t)-1
        endif else begin
         if (i eq max_sec) then begin
;#; Marco, dicembre 2018: metto mint = 1 perché all'inizio di ogni file specyl_?_all.dat scrivo l'ultimo istante del file specyl_?-1_all.dat, e quindi è un duplicato da saltare, specie se uso ibprofiles per calcolare quantità sensibili come la funzione di flusso elicoidale
          mint = 1
          if (n_elements(t) lt end_time) then begin
           maxt = n_elements(t)-1
          endif else begin
           maxt = end_time
          endelse
         endif else begin
          mint = 1
;#; Marco, added
          maxt = n_elements(t)-1
          print,'prova maxt',maxt
         endelse
;#; Marco, marzo 2018, 
         maxt = n_elements(t)-1
        endelse
        print,file
        print,'nt=',nt
        help,bz

        print,'mint ',mint,maxt
;        rpbz = reform(2.*real_part(bz[0:nt[i-1]-1,*,jmin:jmax]))
;        rpbt = reform(2.*real_part(bt[0:nt[i-1]-1,*,jmin:jmax]))
;        rpbr = reform(2.*real_part(br[0:nt[i-1]-1,*,jmin:jmax]))
;        ipbz = reform(2.*imaginary(bz[0:nt[i-1]-1,*,jmin:jmax]))
;        ipbt = reform(2.*imaginary(bt[0:nt[i-1]-1,*,jmin:jmax]))
;        ipbr = reform(2.*imaginary(br[0:nt[i-1]-1,*,jmin:jmax]))
        pbz = reform(2.*(bz[mint:maxt,*,jmin:jmax]))
        pbt = reform(2.*(bt[mint:maxt,*,jmin:jmax]))
        pbr = reform(2.*(br[mint:maxt,*,jmin:jmax]))
        pvz = reform(2.*(vz[mint:maxt,*,jmin:jmax]))
        pvt = reform(2.*(vt[mint:maxt,*,jmin:jmax]))
        pvr = reform(2.*(vr[mint:maxt,*,jmin:jmax]))
        pt = reform(t(mint:maxt))
        help,pt
;        rpvz = reform(2.*real_part(vz[0:nt[i-1]-1,*,jmin:jmax]))
;        rpvt = reform(2.*real_part(vt[0:nt[i-1]-1,*,jmin:jmax]))
;        rpvr = reform(2.*real_part(vr[0:nt[i-1]-1,*,jmin:jmax]))
;        ipvz = reform(2.*imaginary(vz[0:nt[i-1]-1,*,jmin:jmax]))
;        ipvt = reform(2.*imaginary(vt[0:nt[i-1]-1,*,jmin:jmax]))
;        ipvr = reform(2.*imaginary(vr[0:nt[i-1]-1,*,jmin:jmax]))
;        help,rpbz
;        help,rpbt
;        help,rpbr

;#; I remove the 2. factor from the equilibrium harmonic
        if (jmin eq 0) then begin
            pbz[*,*,0] = 0.5*pbz[*,*,0]
            pbt[*,*,0] = 0.5*pbt[*,*,0]
            pbr[*,*,0] = 0.5*pbr[*,*,0]
            pvz[*,*,0] = 0.5*pvz[*,*,0]
            pvt[*,*,0] = 0.5*pvt[*,*,0]
            pvr[*,*,0] = 0.5*pvr[*,*,0]
        endif

        if (i eq min_sec) then begin
            pdbz = temporary(pbz)
            pdbt = temporary(pbt)
            pdbr = temporary(pbr)
            pdvz = temporary(pvz)
            pdvt = temporary(pvt)
            pdvr = temporary(pvr)
            pdt = temporary(pt)

;            print,'memory1',memory(/current)
        endif else begin

            pdbz=[pdbz,temporary(pbz[0:*,*,*])]
            pdbt=[pdbt,temporary(pbt[0:*,*,*])]
            pdbr=[pdbr,temporary(pbr[0:*,*,*])]
            pdvz=[pdvz,temporary(pvz[0:*,*,*])]
            pdvt=[pdvt,temporary(pvt[0:*,*,*])]
            pdvr=[pdvr,temporary(pvr[0:*,*,*])]
            pdt = [pdt,temporary(pt[0:*])]
;            print,'memory2',memory(/current)

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
        
;            print,'memory3',memory(/current)
;    ntp = n_elements(prot)
;    print,'ntp: ',ntp
    print,min(pdt),max(pdt)

;#; Marco, dicembre 2018: tolto per avere sempre ibprofiles.sav aggiornato
;   if (mf ne 1) then begin
   print,'Sto salvando il file.. attendere.. imodeprofiles'
   print,'Attento che c è già il fattore 2 antifft'
   save_file = path+'dat/ibprofiles.sav'
   save,pdt,pror1,pror2,pdbz,pdbt,pdbr,filename=save_file
   save_file = path+'dat/ivprofiles.sav'
   save,pdt,pror1,pror2,pdvz,pdvt,pdvr,filename=save_file
;   endif
   
;#; save the field in modulus&phase notation
   if (mf eq 1) then begin
    bfield_mf_out=path+'dat/imf_bprofiles.sav'
    call_procedure,'mf1',pdt,pdbr,pdbt,pdbz,bfield_mf_out

    vfield_mf_out=path+'dat/imf_vprofiles.sav'
    call_procedure,'mf1',pdt,pdvr,pdvt,pdvz,vfield_mf_out
   endif
   
;#; save the edge field
;   redgebr=reform(pdbr[*,95:100,*])
;   redgebt=reform(pdbt[*,95:100,*])
;   redgebz=reform(pdbz[*,95:100,*])
;   save,filename=path+'dat/ibfield_edge.sav',prot,pror1,pror2,redgebr,redgebt,redgebz

end
