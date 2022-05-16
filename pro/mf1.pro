pro mf1,t, pr,pt,pz,fileout

;#; write the input field using the modulus-phase convention
   mf_z=make_array(n_elements(real_part(pz[*,0,0])),n_elements(real_part(pz[0,*,0])),n_elements(real_part(pz[0,0,*])),/double,value=0.,2)
   mf_t=make_array(n_elements(real_part(pt[*,0,0])),n_elements(real_part(pt[0,*,0])),n_elements(real_part(pt[0,0,*])),/double,value=0.,2)
   mf_r=make_array(n_elements(real_part(pr[*,0,0])),n_elements(real_part(pr[0,*,0])),n_elements(real_part(pr[0,0,*])),/double,value=0.,2)


   mf_z[*,*,*,0] = sqrt( real_part(pz[*,*,*])^2.d + imaginary(pz[*,*,*])^2.d )
   mf_t[*,*,*,0] = sqrt( real_part(pt[*,*,*])^2.d + imaginary(pt[*,*,*])^2.d )
   mf_r[*,*,*,0] = sqrt( real_part(pr[*,*,*])^2.d + imaginary(pr[*,*,*])^2.d )

   if (n_elements(real_part(pr[0,*,0])) eq 101 ) then begin
    pror1 = dindgen(n_elements(real_part(pr[0,*,0]))) / (n_elements(real_part(pr[0,*,0])) - 1)
   endif
   if (n_elements(real_part(pr[0,*,0])) eq 102 ) then begin
    pror2 = dindgen(n_elements(real_part(pr[0,*,0]))) / (n_elements(real_part(pr[0,*,0])) - 2) - 1.d / (2.* (n_elements(real_part(pr[0,*,0])) - 2) )
   endif
   if (n_elements(real_part(pt[0,*,0])) eq 101 ) then begin
    pror1 = dindgen(n_elements(real_part(pt[0,*,0]))) / (n_elements(real_part(pt[0,*,0])) - 1)
   endif
   if (n_elements(real_part(pt[0,*,0])) eq 102 ) then begin
    pror2 = dindgen(n_elements(real_part(pt[0,*,0]))) / (n_elements(real_part(pt[0,*,0])) - 2) - 1.d / (2.* (n_elements(real_part(pt[0,*,0])) - 2) )
   endif
   
   mf_z[*,*,*,1] = atan( imaginary(pz[*,*,*]),real_part(pz[*,*,*]) )
  ;#; atan correction to have phases between 0,2*pi
   aa = where(reform(mf_z[*,*,*,1]) lt 0.)
   if (total(aa) gt 0.) then begin
    aaa=array_indices(reform(mf_z[*,*,*,1]),aa)
    help,aaa
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
   save,filename=fileout,t,mf_r,mf_t,mf_z,pror1,pror2

end

