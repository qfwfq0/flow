pro isavelprofiles,path,mf,d3,tmin,tmax
   path0 = path
   RR = read_major_radius(path)

;   restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
   restore_ifile = path+'dat/ivprofiles.sav'
   ex_file = file_test(restore_ifile)
   if (ex_file ne 1) then begin
    print,'need to create the fields file ivprofiles.sav'
    call_procedure,'isavemodeprofiles',path,'4',d3
   endif
   restore,restore_ifile
   
   ntp = fix(tmax) - fix(tmin)
   jmax = n_elements(pdvr[0,0,*])-1
   jmax = 2
   nx = n_elements(pdvr[0,*,0])-1   
   print,'nx',nx
   help,pdvt,pdvz
   dr = 1.d / nx
    
;#; declaration of the components of angular momentum
   llr=make_array(ntp,jmax+1,nx+2,/dcomplex,value=0.)
   llt=make_array(ntp,jmax+1,nx+2,/dcomplex,value=0.)
   llz=make_array(ntp,jmax+1,nx+2,/dcomplex,value=0.)

   print,'I need to compute the angular momentum for '+strtrim(string(ntp,'(i0)'))+' time steps, from t='+strtrim(tmin)+' until t='+strtrim(tmax)

   
   for it=fix(tmin),fix(tmax)-1 do begin
     if (it mod 20 eq 0) then begin
      print,'computing angular momentum for itp=',it
     endif
       for j=0,jmax do begin
        mn = mnum(j,nz,nanf,d3)
;#; first interpolation of vr on the x2 mesh
        dvr = interpol(reform(pdvr[it,*,j]), pror1, pror2 )
        if (mn[0] eq 0) then begin
         dvr[0] = -dvr[1]
        endif else begin
         if (mn[1] eq 1) then begin
          dvr[0] = dvr[1]
         endif else begin
          dvr[0] = -dvr[1]
         endelse
        endelse

        for i=1,nx+1 do begin
          llr[it,j,i] = dcomplex(0.,0.)
          llt[it,j,i] = dcomplex(- pror2[i]*real_part(pdvz[it,i,j]) , - pror2[i]*imaginary(pdvz[it,i,j]))
          llz[it,j,i] = dcomplex(pror2[i]*real_part(pdvt[it,i,j]), -pror2[i]*imaginary(pdvt[it,i,j]))
       endfor
;;#; regularity conditions at the core
           if (mn[0] eq 0) then begin
            llz[it,j,0] = llz[it,j,1]
            llt[it,j,0] = -llt[it,j,1]
            llr[it,j,0] = -llr[it,j,1]
           endif else begin
            if (mn[0] eq 1) then begin
             llz[it,j,0] = -llz[it,j,1]
             llt[it,j,0] = +llt[it,j,1]
             llr[it,j,0] = +llr[it,j,1]
            endif else begin
             llz[it,j,0] = -llz[it,j,1]
             llt[it,j,0] = -llt[it,j,1]
             llr[it,j,0] = -llr[it,j,1]
            endelse
           endelse
       endfor
   endfor

   save_file = path+'dat/ilprofiles.sav'
   save,filename=save_file,pror2,pdt,llr,llt,llz,nx

;#; save the field in modulus&phase notation
;   if (mf eq 1) then begin
;    jfield_mf_out=path+'dat/imf_lprofiles.sav'
;    call_procedure,'mf1',pdt,pjr,pjt,pjz,jfield_mf_out
;   endif
end
