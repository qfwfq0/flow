pro energy,t,br,bt,bz,str,nsec,path

;#; define major radius
;RR=read_major_radius(path)
RR=4.d
 print,'RR_energia',RR
fact = 4. * !Dpi * !Dpi  * RR[0]

print,fact
;#; define radial profiles
lx = n_elements(br[0,*,0]) * 1.d
pror1 = dindgen(lx) / (lx-1) 
pror2 = dindgen(lx + 1) / (lx-1) + 1.d / (2.d*(lx-1))

comp = make_array(n_elements(br[*,0,0]),n_elements(br[0,0,*]),/double,value=0.)
total_en = make_array(n_elements(br[*,0,0]),/double,value=0.)

;#; understand which kind of energy is being computed'
case str of
 'b': mssg='Computing magnetic energy at time'
 'v': mssg='Computing kynetic energy at time'
endcase

for tq=0, n_elements(br[*,0,0])-1 do begin

if (tq mod 50) eq 0 then print, mssg,t(tq)
for j=0, n_elements(br[0,0,*])-1 do begin
  if (j eq 0) then begin
   comp(tq,j) = int_tabulated(pror1,pror1*real_part(br[tq,*,j])^2.) + int_tabulated(pror2,pror2*real_part(bt[tq,*,j])^2.) + int_tabulated(pror2,pror2*real_part(bz[tq,*,j])^2.)
   total_en[tq] = fact * comp(tq,j)
  endif else begin
   comp(tq,j) = 2. * int_tabulated(pror1,pror1*real_part(br[tq,*,j])^2.) + 2. * int_tabulated(pror1,pror1*imaginary(br[tq,*,j])^2.) + 2. * int_tabulated(pror2,pror2*real_part(bt[tq,*,j])^2.) + 2. * int_tabulated(pror2,pror2*imaginary(bt[tq,*,j])^2.) + 2.*int_tabulated(pror2,pror2*real_part(bz[tq,*,j])^2.) + 2.*int_tabulated(pror2,pror2*imaginary(bz[tq,*,j])^2.)
   total_en[tq] = total_en[tq] + fact * comp(tq,j)
  endelse
;#; end of the cycle in j (janz)  
endfor
   total_en[tq] = 0.5 * total_en[tq]
;#; end of the cycle in tq (time)  
endfor

;#; add the fact also to the comp vector
;#; Marco, 20/03/2017: after benchmarking with old specyl-version: I did not put a 0.5 factor in the comp vector
 comp = 0.5 * comp * fact
filename = path+'dat/'+str+'_energy_'+strtrim(string(nsec,format='(i0)'))+'.sav'
print,filename
save,filename=filename,t,comp,total_en

end

