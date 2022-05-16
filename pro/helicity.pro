pro helicity,path

;#; define major radius
RR=read_major_radius(path)
;RR=4.d
 print,'RR_energia',RR
fact = 4. * !Dpi * !Dpi  * RR[0]
print,fact

restore,path+'dat/ibprofiles.sav'
restore,path+'dat/iaprofiles.sav'
;#; define radial profiles
lx = n_elements(pdbr[0,*,0]) * 1.d
pror1 = dindgen(lx) / (lx-1) 
pror2 = dindgen(lx + 1) / (lx-1) + 1.d / (2.d*(lx-1))

a_comp = make_array(n_elements(pdbr[*,0,0]),n_elements(pdbr[0,0,*]),/double,value=0.)
total_hel = make_array(n_elements(pdbr[*,0,0]),/double,value=0.)

;#; understand which kind of energy is being computed'
mssg='Computing magnetic helicity at time'

for tq=0, n_elements(pdbr[*,0,0])-1 do begin

if (tq mod 50) eq 0 then print, mssg,pdt(tq)
for j=0, n_elements(pdbr[0,0,*])-1 do begin
  if (j eq 0) then begin
   a_comp(tq,j) = int_tabulated(pror2,pror2*real_part(pdbt[tq,*,j])*real_part(pdat[tq,*,j])) + int_tabulated(pror2,pror2*real_part(pdbz[tq,*,j])*real_part(pdaz[tq,*,j]))
   total_hel[tq] = fact * a_comp(tq,j)
  endif else begin
   a_comp(tq,j) = 2. * int_tabulated(pror2,pror2*real_part(pdbt[tq,*,j])*real_part(pdat[tq,*,j])) + 2. * int_tabulated(pror2,pror2*imaginary(pdbt[tq,*,j])*imaginary(pdat[tq,*,j])) + 2.*int_tabulated(pror2,pror2*real_part(pdbz[tq,*,j])*real_part(pdaz[tq,*,j])) + 2.*int_tabulated(pror2,pror2*imaginary(pdbz[tq,*,j])*imaginary(pdaz[tq,*,j]))
   total_hel[tq] = total_hel[tq] + fact * a_comp(tq,j)
  endelse
;#; end of the cycle in j (janz)  
endfor
   total_hel[tq] = 0.5 * total_hel[tq]
;#; end of the cycle in tq (time)  
endfor

;#; add the fact also to the comp vector
;#; Marco, 20/03/2017: after benchmarking with old specyl-version: I did not put a 0.5 factor in the comp vector
 a_comp = 0.5 * a_comp * fact
 filename = path+'dat/helicity.sav'
 print,filename
 save,filename=filename,tt,a_comp,total_hel

end

