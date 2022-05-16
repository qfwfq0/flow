pro isavetauprofiles,path,mf,d3,tmin,tmax
;#; this pro computes the momentum of various forces acting on the plasma (JxB, nu*lapl(v) etc etc)
   path0 = path
   RR = read_major_radius(path)

;   restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
   restore_ifile = path+'dat/ijprofiles.sav'
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
    
;#; JxB
;#; declaration of the components of momentum
   jbtaur=make_array(ntp,jmax+1,nx+2,/dcomplex,value=0.)
   jbtaut=make_array(ntp,jmax+1,nx+2,/dcomplex,value=0.)
   jatauz=make_array(ntp,jmax+1,nx+2,/dcomplex,value=0.)

   print,'I need to compute the JxB momentum for '+strtrim(string(ntp,'(i0)'))+' time steps, from t='+strtrim(tmin)+' until t='+strtrim(tmax)

   
   for it=fix(tmin),fix(tmax)-1 do begin
     if (it mod 20 eq 0) then begin
      print,'computing angular momentum for itp=',it
     endif
       for j=0,jmax do begin
        mn = mnum(j,nz,nanf,d3)
;#; interpolation of jt,jz,br on the x2 mesh
        djt = interpol(reform(pjt[it,*,j]), pror1, pror2 )
        djz = interpol(reform(pjt[it,*,j]), pror1, pror2 )
        dbr = interpol(reform(pdbr[it,*,j]), pror1, pror2 )
        if (mn[0] eq 0) then begin
         dbr[0] = - dbr[1]
         djt[0] = - djt[1]
         djz[0] =   djz[1]
        endif else begin
         if (mn[1] eq 1) then begin
          dbr[0] =  dbr[1]
          djt[0] =  djt[1]
          djz[0] = -djz[1]
         endif else begin
          dbr[0] = - dbr[1]
          djt[0] = - djt[1]
          djz[0] = - djz[1]
         endelse
        endelse

        for i=1,nx+1 do begin
          jbtaur[it,j,i] = dcomplex(0.,0.)
          jbtaut[it,j,i] = dcomplex( pror2[i]*real_part(pjr[it,i,j])*real_part() , - pror2[i]*imaginary(pdvz[it,i,j]))
          jbtauz[it,j,i] = dcomplex(pror2[i]*real_part(pdvt[it,i,j]), -pror2[i]*imaginary(pdvt[it,i,j]))
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
