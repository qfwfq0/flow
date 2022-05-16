pro mf_s,t, pr,pt,pz,fileout,str

;#; write the input field using the modulus-phase convention
   mf_z=make_array(n_elements(real_part(pz[*,0,0])),n_elements(real_part(pz[0,*,0])),n_elements(real_part(pz[0,0,*])),/double,value=0.,2)
   mf_t=make_array(n_elements(real_part(pt[*,0,0])),n_elements(real_part(pt[0,*,0])),n_elements(real_part(pt[0,0,*])),/double,value=0.,2)
   mf_r=make_array(n_elements(real_part(pr[*,0,0])),n_elements(real_part(pr[0,*,0])),n_elements(real_part(pr[0,0,*])),/double,value=0.,2)


   mf_z[*,*,*,0] = sqrt( real_part(pz[*,*,*])^2.d + imaginary(pz[*,*,*])^2.d )
   mf_t[*,*,*,0] = sqrt( real_part(pt[*,*,*])^2.d + imaginary(pt[*,*,*])^2.d )
   mf_r[*,*,*,0] = sqrt( real_part(pr[*,*,*])^2.d + imaginary(pr[*,*,*])^2.d )


   mf_z[*,*,*,1] = atan( imaginary(pz[*,*,*]),real_part(pz[*,*,*]) )
  ;#; atan correction to have phases between 0,2*pi
   aa = where(reform(mf_z[*,*,*,1]) lt 0.)
   if (total(aa) gt 0.) then begin
    aaa=array_indices(reform(mf_z[*,*,*,1]),aa)
    for k=0,n_elements(aaa[0,*])-1 do begin
     mf_z[aaa[0,k],aaa[1,k],aaa[2,k],1] = mf_z[aaa[0,k],aaa[1,k],aaa[2,k],1] + 2.*!dpi
    endfor
   endif
  
   mf_t[*,*,*,1] = atan( imaginary(pt[*,*,*]),real_part(pt[*,*,*]) )
  ;#; atan correction to have phases between 0,2*pi
   aa = where(reform(mf_t[*,*,*,1]) lt 0.)
   if (total(aa) gt 0.) then begin
    aaa=array_indices(reform(mf_t[*,*,*,1]),aa)
    for k=0,n_elements(aaa[0,*])-1 do begin
     mf_t[aaa[0,k],aaa[1,k],aaa[2,k],1] = mf_t[aaa[0,k],aaa[1,k],aaa[2,k],1] + 2.*!dpi
    endfor
   endif
  
   mf_r[*,*,*,1] = atan( imaginary(pr[*,*,*]), real_part(pr[*,*,*]) )
  ;#; atan correction to have phases between 0,2*pi
   aa = where(reform(mf_r[*,*,*,1]) lt 0.)
   if (total(aa) gt 0.) then begin
    aaa=array_indices(reform(mf_r[*,*,*,1]),aa)
    for k=0,n_elements(aaa[0,*])-1 do begin
     mf_r[aaa[0,k],aaa[1,k],aaa[2,k],1] = mf_r[aaa[0,k],aaa[1,k],aaa[2,k],1] + 2.*!dpi
    endfor
   endif

   print,'saving...'
   if (str ne '') then begin
    if (str eq 'el') then begin
     mf_eq_elr = reform(mf_r) & mf_eq_elt = reform(mf_t) & mf_eq_elz = reform(mf_z)
     save,filename=fileout,t,mf_eq_elr,mf_eq_elt,mf_eq_elz
    endif
   endif else begin
    print,'saving...'
    save,filename=fileout,t,mf_r,mf_t,mf_z
   endelse


end

