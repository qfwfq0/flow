function mf, pr,pt,pz
;#; write the input field using the modulus-phase convention
 mf_z=make_array(n_elements(pz[*,0,0]),n_elements(rpz[0,*,0]),n_elements(rpz[0,0,*]),/double,value=0.)
 mf_t=make_array(2,n_elements(rpt[*,0,0]),n_elements(rpt[0,*,0]),n_elements(rpt[0,0,*]),/double,value=0.)
 mf_r=make_array(2,n_elements(rpr[*,0,0]),n_elements(rpr[0,*,0]),n_elements(rpr[0,0,*]),/double,value=0.)
 mf_z[0,*,*,*] = sqrt( rpz[*,*,*]^2.d + ipz[*,*,*]^2.d )
 mf_t[0,*,*,*] = sqrt( rpt[*,*,*]^2.d + ipt[*,*,*]^2.d )
 mf_r[0,*,*,*] = sqrt( rpr[*,*,*]^2.d + ipr[*,*,*]^2.d )
 mf_z[1,*,*,*] = atan( ipz[*,*,*],rpz[*,*,*] )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_z[1,*,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_z[1,*,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_z[1,aaa[0,k],aaa[1,k],aaa[2,k]] = mf_z[1,aaa[0,k],aaa[1,k],aaa[2,k]] + 2.*!dpi
  endfor
 endif

 mf_t[1,*,*,*] = atan( rpt[*,*,*],ipt[*,*,*] )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_t[1,*,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_t[1,*,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_t[1,aaa[0,k],aaa[1,k],aaa[2,k]] = mf_t[1,aaa[0,k],aaa[1,k],aaa[2,k]] + 2.*!dpi
  endfor
 endif

 mf_r[1,*,*,*] = atan( rpr[*,*,*],ipr[*,*,*] )
;#; atan correction to have phases between 0,2*pi
 aa = where(reform(mf_r[1,*,*,*]) lt 0.)
 if (total(aa) gt 0.) then begin
  aaa=array_indices(reform(mf_r[1,*,*,*]),aa)
  for k=0,n_elements(aaa[0,*])-1 do begin
   mf_r[1,aaa[0,k],aaa[1,k],aaa[2,k]] = mf_r[1,aaa[0,k],aaa[1,k],aaa[2,k]] + 2.*!dpi
  endfor
 endif
 
 save,filename=fileout,dt,mf_r,mf_t,mf_z

end

