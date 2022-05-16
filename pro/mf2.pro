pro mf2,t, pr,pt,pz,fileout
;#! VERY IMPORTANT: this version keeps into account that SpeCyl may write, as an output, 
;SpeCyl che vuole il campo reale corrispondente scritto come
; 2.* Re(x) * cos(u) - 2. * Im(x) * sin(u)
;#; write the input field using the modulus-phase convention
   mf_z=make_array(n_elements(real_part(pz[*,0,0])),n_elements(real_part(pz[0,*,0])),n_elements(real_part(pz[0,0,*])),/double,value=0.,2)
   mf_t=make_array(n_elements(real_part(pt[*,0,0])),n_elements(real_part(pt[0,*,0])),n_elements(real_part(pt[0,0,*])),/double,value=0.,2)
   mf_r=make_array(n_elements(real_part(pr[*,0,0])),n_elements(real_part(pr[0,*,0])),n_elements(real_part(pr[0,0,*])),/double,value=0.,2)


   mf_z[*,*,*,0] = sqrt( real_part(pz[*,*,*])^2.d + imaginary(pz[*,*,*])^2.d )
   mf_t[*,*,*,0] = sqrt( real_part(pt[*,*,*])^2.d + imaginary(pt[*,*,*])^2.d )
   mf_r[*,*,*,0] = sqrt( real_part(pr[*,*,*])^2.d + imaginary(pr[*,*,*])^2.d )


;   mf_z[*,*,*,1] = atan( imaginary(pz[*,*,*]),real_part(pz[*,*,*]) )
   mf_z[*,*,*,1] = atan(-imaginary(pz[*,*,*]),real_part(pz[*,*,*]) )
  ;#; atan correction to have phases between 0,2*pi
   aa = where(reform(mf_z[*,*,*,1]) lt 0.)
   if (total(aa) gt 0.) then begin
    aaa=array_indices(reform(mf_z[*,*,*,1]),aa)
    for k=0,n_elements(aaa[0,*])-1 do begin
     mf_z[aaa[0,k],aaa[1,k],aaa[2,k],1] = mf_z[aaa[0,k],aaa[1,k],aaa[2,k],1] + 2.*!dpi
    endfor
   endif
  
;   mf_t[*,*,*,1] = atan( imaginary(pt[*,*,*]),real_part(pt[*,*,*]) )
   mf_t[*,*,*,1] = atan( -imaginary(pt[*,*,*]),real_part(pt[*,*,*]) )
;   print,pt[0,101,1]
;   print, 'a',atan( -imaginary(pt[0,101,1]),real_part(pt[0,101,1]) )
;   print, 'b',atan( imaginary(pt[0,101,1]),real_part(pt[0,101,1]) )
;   print,mf_t[0,101,1,1]
  ;#; atan correction to have phases between 0,2*pi
   aa = where(reform(mf_t[*,*,*,1]) lt 0.)
   if (total(aa) gt 0.) then begin
    aaa=array_indices(reform(mf_t[*,*,*,1]),aa)
    for k=0,n_elements(aaa[0,*])-1 do begin
     mf_t[aaa[0,k],aaa[1,k],aaa[2,k],1] = mf_t[aaa[0,k],aaa[1,k],aaa[2,k],1] + 2.*!dpi
    endfor
   endif
;   print,mf_t[0,101,1,1]
  
;   mf_r[*,*,*,1] = atan( imaginary(pr[*,*,*]), real_part(pr[*,*,*]) )
   mf_r[*,*,*,1] = atan( -imaginary(pr[*,*,*]), real_part(pr[*,*,*]) )
  ;#; atan correction to have phases between 0,2*pi
   aa = where(reform(mf_r[*,*,*,1]) lt 0.)
   if (total(aa) gt 0.) then begin
    aaa=array_indices(reform(mf_r[*,*,*,1]),aa)
    for k=0,n_elements(aaa[0,*])-1 do begin
     mf_r[aaa[0,k],aaa[1,k],aaa[2,k],1] = mf_r[aaa[0,k],aaa[1,k],aaa[2,k],1] + 2.*!dpi
    endfor
   endif
   
   print,'saving...'
   save,filename=fileout,t,mf_r,mf_t,mf_z

end

