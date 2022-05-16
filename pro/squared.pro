pro squared,path,vect,d3

;path='./'
;vect='b'
;d3=1
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse

if (vect eq 'b') then begin
;#; magnetic field
   restore_ifile = path+'dat/ibprofiles.sav'
   ex_file = file_test(restore_ifile)
   if (ex_file ne 1) then begin
    print,'need to create the fields file ibprofiles.sav'
    call_procedure,'isavemodeprofiles',path,'0',d3
   endif
   restore,restore_ifile

endif else begin
;#; velocity field
   restore_ifile = path+'dat/ivprofiles.sav'
   ex_file = file_test(restore_ifile)
   if (ex_file ne 1) then begin
    print,'need to create the fields file ivprofiles.sav'
    call_procedure,'isavemodeprofiles',path,'0',d3
   endif
   restore,restore_ifile
endelse

if (vect eq 'b') then begin

 square = pdbt*dcomplex(0.d,0.d)
 help,square
 for it=0, n_elements(square[*,0,0]) -1 do begin
  if (it mod 50 eq 0) then print, 'computing itp=',it
  for j=0, n_elements(square[0,0,*]) -1  do begin
   ; define mn! 
     mn = mnum(j,nz,nanf,d3)
;     print,'squared',mn
     cnvl = convo(path,mn[0],mn[1],d3)
    for p = 0, n_elements(cnvl[0,*])-1 do begin
     j1 = cnvl[0,p]
     j2 = cnvl[1,p]
;#; first interpolation of the required br on the x2 mesh
     dbr1 = interpol(reform(pdbr[it,*,j1]), pror1, pror2 )
     dbr2 = interpol(reform(pdbr[it,*,j2]), pror1, pror2 )
     mn1 = mnum(j1,nz,nanf,d3)
     mn2 = mnum(j2,nz,nanf,d3)
     if (mn1[0] eq 0) then begin
      dbr1[0] = -dbr1[1]
     endif else begin
      if (mn1[1] eq 1) then begin
       dbr1[0] = dbr1[1]
      endif else begin
       dbr1[0] = -dbr1[1]
      endelse
     endelse
     if (mn2[0] eq 0) then begin
      dbr2[0] = -dbr2[1]
     endif else begin
      if (mn2[1] eq 1) then begin
       dbr2[0] = dbr2[1]
      endif else begin
       dbr2[0] = -dbr2[1]
      endelse
     endelse

;#; definition of the multiplication factor
     fact = mlpf(mn1,mn2)
;     print,'j1',j1,'j2',j2,fact

     square(it,*,j) = square(it,*,j) + fact*dbr1*dbr2 + fact*pdbt[it,*,j1]*pdbt[it,*,j2] + fact*pdbz[it,*,j1]*pdbz[it,*,j2]
;#; Marco, end of the convolution cycle
    endfor
;#; Marco, end of the cycle on mode number
   endfor
;#; Marco, end of the cycle on time
  endfor
   save_ifile = path+'dat/ibsq_profiles.sav'
   save,filename=save_ifile,square,pror2,pdt
;#; Marco, end of the 'b' case, begins 'v' case
endif else begin

 square = pdvt*dcomplex(0.d,0.d)
 help,square
 for it=0, n_elements(square[*,0,0]) -1 do begin
  if (it mod 50 eq 0) then print, 'computing itp=',it
  for j=0, n_elements(square[0,0,*]) -1  do begin
   ; define mn! 
     mn = mnum(j,nz,nanf,d3)
     cnvl = convo(path,mn[0],mn[1],d3)
    for p = 0, n_elements(cnvl[0,*])-1 do begin
     j1 = cnvl[0,p]
     j2 = cnvl[1,p]
;#; first interpolation of the required vr on the x2 mesh
     dbv1 = interpol(reform(pdvr[it,*,j1]), pror1, pror2 )
     dbv2 = interpol(reform(pdvr[it,*,j2]), pror1, pror2 )
     mn1 = mnum(j1,nz,nanf,d3)
     mn2 = mnum(j2,nz,nanf,d3)
     if (mn1[0] eq 0) then begin
      dvr1[0] = -dvr1[1]
     endif else begin
      if (mn1[1] eq 1) then begin
       dvr1[0] = dvr1[1]
      endif else begin
       dvr1[0] = -dvr1[1]
      endelse
     endelse
     if (mn2[0] eq 0) then begin
      dvr2[0] = -dvr2[1]
     endif else begin
      if (mn2[1] eq 1) then begin
       dvr2[0] = dvr2[1]
      endif else begin
       dvr2[0] = -dvr2[1]
      endelse
     endelse

;#; definition of the multiplication factor
    fact = mlpf(mn1,mn2)

     square(it,*,j) = square(it,*,j) + fact*dvr1*dvr2 + fact*pdvt[it,*,j1]*pdvt[it,*,j2] + fact*pdvz[it,*,j1]*pdvz[it,*,j2]
;#; Marco, end of the convolution cycle
    endfor
;#; Marco, end of the cycle on mode number
   endfor
;#; Marco, end of the cycle on time
  endfor
   save_ifile = path+'dat/ivsq_profiles.sav'
   save,filename=save_ifile,square,pror2,pdt
endelse
end
