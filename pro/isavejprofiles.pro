pro isavejprofiles,path,mf,d3
   path0 = path
   RR = read_major_radius(path)

   print,'ciao'
   result = file_test(path+'dat/spectrum.sav')
   if (result ne 1) then begin
    call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
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
   
   ntp = n_elements(pdt)
   jmax = n_elements(pdbr[0,0,*])-1
   if (d3 eq 1) then begin
    jmax=n_elements(pdbr[0,0,*])-1
   endif else begin
    jmax=5
   endelse
   print,'jmax= ',jmax
   lx = n_elements(pdbr[0,*,0])-1   
   dr = 1.d / lx
    
;;#; declaration of real part and imaginary part of the components of current density
;   ipjr = imaginary(pdbz)*0.d
;   rpjr = real_part(pdbz)*0.d
;   ipjt = imaginary(pdbr)*0.d
;   rpjt = real_part(pdbr)*0.d
;   ipjz = imaginary(pdbr)*0.d
;   rpjz = real_part(pdbr)*0.d

;#; declaration of the components of current density
   pjr=make_array(ntp,n_elements(pdbz[0,*,0]),jmax+1,/dcomplex,value=0.)
   pjt=make_array(ntp,n_elements(pdbr[0,*,0]),jmax+1,/dcomplex,value=0.)
   pjz=make_array(ntp,n_elements(pdbr[0,*,0]),jmax+1,/dcomplex,value=0.)

   print,'I need to compute the current for '+string(ntp,format='(i0)')+' time steps'

   for it=0,ntp-1 do begin
     if (it mod 50 eq 0) then begin
      print,'computing current for itp=',it
     endif
       for j=0,jmax do begin
;           mn = mnum(j,nz,nanf,d3)
;#; Marco, 31/10/2018: ho avuto problemi con j=135 (cioè l'ultimo elemento dello spettro).
;#; sostituisco nz con nzcon che ha un elemento in più
           mn = mnum(j,nzcon,nanf,d3)
;           if (d3 eq 0) then print,'mn 2d',mn
           for ix=0,lx do begin
               pjz[it,ix,j] = dcomplex( 1./(pror1(ix)*dr)*(pror2(ix+1)*reform(real_part(pdbt[it,ix+1,j]))-pror2(ix)*reform(real_part(pdbt[it,ix,j]))) + mn(0)/pror1(ix)*reform(imaginary(pdbr[it,ix,j])) ,  1./(pror1(ix)*dr)*(pror2(ix+1)*reform(imaginary(pdbt[it,ix+1,j]))-pror2(ix)*reform(imaginary(pdbt[it,ix,j]))) - mn(0)/pror1(ix)*reform(real_part(pdbr[it,ix,j])) )
               pjt[it,ix,j] = dcomplex( -1./dr*reform(real_part(pdbz[it,ix+1,j])-real_part(pdbz[it,ix,j])) - mn(1)/RR(0)*reform(imaginary(pdbr[it,ix,j])) , -1./dr*reform(imaginary(pdbz[it,ix+1,j])-imaginary(pdbz[it,ix,j])) + mn(1)/RR(0)*reform(real_part(pdbr[it,ix,j])) )
           endfor
           for ix=0,lx+1 do begin
               pjr[it,ix,j] = dcomplex ( -mn(0)/pror2(ix)*reform(imaginary(pdbz[it,ix,j])) + mn(1)/RR(0)*reform(imaginary(pdbt[it,ix,j])) ,  +mn(0)/pror2(ix)*reform(real_part(pdbz[it,ix,j])) - mn(1)/RR(0)*reform(real_part(pdbt[it,ix,j])) )
           endfor
;#; regularity conditions at the core
           if (mn[0] eq 0) then begin
            pjz[it,0,j] = pjz[it,1,j]
            pjt[it,0,j] = 0.d
            pjr[it,0,j] = -pjr[it,1,j]
           endif else begin
            if (mn[0] eq 1) then begin
             pjz[it,0,j] = 0.d
;             pjt[it,0,j] = pjt[it,1,j]
;             pjr[it,0,j] = pjr[it,1,j]
            endif else begin
             pjz[it,0,j] = 0.d
             pjt[it,0,j] = 0.d
             pjr[it,0,j] = -pjr[it,1,j]
            endelse
           endelse
       endfor
   endfor

   save_file = path+'dat/ijprofiles.sav'
   save,filename=save_file,pror1,pror2,pdt,pjr,pjt,pjz

;#; save the field in modulus&phase notation
   if (mf eq 1) then begin
    jfield_mf_out=path+'dat/imf_jprofiles.sav'
    call_procedure,'mf1',pdt,pjr,pjt,pjz,jfield_mf_out
   endif
end
