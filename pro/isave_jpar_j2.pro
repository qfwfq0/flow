pro isave_jpar_j2, path, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax

   mf = 0
   path0 = path
   RR = read_major_radius(path)
   dissipation = read_dissipation(path)
   print,'dissipation',dissipation.eta,dissipation.nu
   if (dissipation.eta*dissipation.nu le 0.d) then begin 
    print,'si è letta male la dissipazione'
    stop
   endif

   result = file_test(path+'dat/spectrum.sav')
;   if (result ne 1) then begin
;    call_procedure, 'read_spectrum',path,my,mm,nanf,nz
;   endif else begin
;    restore,path+'dat/spectrum.sav'
;   endelse
   if (result ne 1) then begin
    call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
    print,'nanf[1]',nanf[1]
   endif else begin
    restore,path+'dat/spectrum.sav'
    nzcon=nz
    nzcon[0]=nz[0]+1
   endelse


   restore_ifile = path+'dat/ibprofiles.sav'
   ex_file = file_test(restore_ifile)
   print,'ex_file',ex_file
   if (ex_file ne 1) then begin
    mf0 = mf
    print,'need to create the fields file ibprofiles.sav'
    print,'d3= ',d3
    mf = 0
    call_procedure,'isavemodeprofiles',path,mf,d3
    mf = mf0
   endif
   restore,restore_ifile
   restore_jfile = path+'dat/ijprofiles.sav'
   ex_file = file_test(restore_jfile)
   print,'ex_file',ex_file
   if (ex_file ne 1) then begin
    mf0 = mf
    print,'need to create the fields file ijprofiles.sav'
    print,'d3= ',d3
    mf = 0
    call_procedure,'isavejprofiles',path,mf,d3
    mf = mf0
   endif
   restore,restore_jfile
   ntp = n_elements(pdt)
   jmax = n_elements(pdbr[0,0,*])-1
   print,'jmax= ',jmax
   lx = n_elements(pdbr[0,*,0])-1   
   dr = 1.d / lx
;#; angular meshes
   vecth = dindgen(nth)/(nth-1) * 2. * !Dpi
   if (nzz eq 1) then begin
    vecz = 0.
   endif else begin
    vecz = dindgen(nzz)/(nzz-1) * 2. * !Dpi * RR[0]
   endelse
 
;#; check input times
   if (itp0 eq !null) then begin
   itp0 = 0
   endif
   if (itp1 eq !null) then begin
   itp1 = n_elements(pdt) - 1
   endif
   if (itp0 eq !null) then begin
    print,'define itp0!'
   endif
   if (itp1 eq -1) then begin
   itp1 = n_elements(pdt) - 1
   help,pdt
   endif
   if (itp0 eq -1) then begin
   itp0 = n_elements(pdt) - 1
   endif
   itp0 = fix(itp0,type=3)
   itp1 = fix(itp1,type=3)
   nitp = itp1 - itp0 + 1
   nsaved = ceil(float(nitp)/deltaitp)
   iitp= indgen(nsaved)*deltaitp + itp0 
   print,deltaitp,itp0,itp1
   help,iitp

;#; ricostruzione dei campi magnetico e di corrente   
   q=0
   help, iitp
   for q = 0 , nsaved-1 do begin
   str='j'
   restore_jreconstruction = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
   ex_file = file_test(restore_jreconstruction)
   if (ex_file ne 1) then begin
     mf0 = mf
     print,'need to create the reconstructed fields file'
     print,'d3= ',d3
     mf = 0
     call_procedure,'antifft_reconst_rfp',path,'j', d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
     call_procedure,'antifft_reconst_rfp',path,'b', d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
     mf = mf0
   endif
   endfor

  
   mesh_size = read_mesh(path)
   nbr = make_array(nsaved,mesh_size+2,nth,nzz,/dcomplex,value=0.d)
   nbt = make_array(nsaved,mesh_size+2,nth,nzz,/dcomplex,value=0.d)
   nbz = make_array(nsaved,mesh_size+2,nth,nzz,/dcomplex,value=0.d)
   njr = make_array(nsaved,mesh_size+2,nth,nzz,/dcomplex,value=0.d)
   njt = make_array(nsaved,mesh_size+2,nth,nzz,/dcomplex,value=0.d)
   njz = make_array(nsaved,mesh_size+2,nth,nzz,/dcomplex,value=0.d)

   print,'second time help ITP'
    for q = 0 , nsaved-1 do begin
     str='j'
     print,'q=',q
;     restore_jreconstruction = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
     restore_jreconstruction = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'_pert.sav'
     restore,restore_jreconstruction
     ;#; per prima cosa interpolo il campo jt sulla mesh x2
     nnjt = interpolate_x2_vec(reform(jt),path,d3) 
     nnjz = interpolate_x2_vec(reform(jz),path,d3) 
     njr[q,*,*,*] = jr
     njt[q,*,*,*] = nnjt
     njz[q,*,*,*] = nnjz

     str='b'
;     restore_breconstruction = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
     restore_breconstruction = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'_pert.sav'
     restore,restore_breconstruction
     ;#; per prima cosa interpolo il campo br sulla mesh x2
     nnbr = interpolate_x2_vec(reform(br),path,d3) 
     nbr[q,*,*,*] = nnbr
     nbt[q,*,*,*] = bt
     nbz[q,*,*,*] = bz

    endfor
;#; declaration of the components of current density
    bsq = make_array(nsaved,mesh_size+2,nth,nzz,/double,value=0.d)
    jsq = make_array(nsaved,mesh_size+2,nth,nzz,/double,value=0.d)
    jpar = make_array(nsaved,mesh_size+2,nth,nzz,/double,value=0.d)

    for q = 0 , nsaved-1 do begin
     bsq(q,*,*,*) = nbr(q,*,*,*)^2.d + nbt(q,*,*,*)^2.d + nbz(q,*,*,*)^2.d
     jsq(q,*,*,*) = njr(q,*,*,*)^2.d + njt(q,*,*,*)^2.d + njz(q,*,*,*)^2.d
     jpar(q,*,*,*) = ( njr(q,*,*,*)*nbr(q,*,*,*) + njt(q,*,*,*)*nbt(q,*,*,*) +njz(q,*,*,*)*nbz(q,*,*,*) ) / bsq(q,*,*,*)
     time_file = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/etaj2_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
     absq = reform(bsq[q,*,*,*])
     ajsq = reform(jsq[q,*,*,*])
     ajpar = reform(jpar[q,*,*,*])
     atime=pdt(itp(q))
     save,filename=time_file,absq,ajsq,ajpar,pror2,vecth,vecz,atime

    endfor
    outfile = path+'dat/jpar_j2_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'_itp'+strtrim(string(itp0,format='(i0)'))+strtrim(string(itp1,format='(i0)'))+'.sav'
    print,outfile   
    time=pdt(iitp)
    save,filename=outfile,bsq,jsq,jpar,pror2,vecth,vecz,time
end

