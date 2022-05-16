function integral2,vect,xx
;#; general purpose function to calculate the integral of a function (vect) along the mesh x2 (vect)
;#; x2=findgen(lx+2)/lx-1./(2.*lx)
;#; the result will be stored in the mesh x1(i) = findgen(lx+1)/lx

 integ = make_array(n_elements(xx)-1,n_elements(vect(0,*)),/double,value=0.)
 yy = findgen(n_elements(xx)-1) / (n_elements(xx)-2)
 dx = 1.d / (n_elements(xx)-2)
 
 for q=0, n_elements(vect(0,*))-1 do begin
  integ(0,q) = 0.d
  for i=1,n_elements(yy)-1 do begin
   integ(i,q) = integ(i-1,q) + vect(i,q) * dx
  endfor
 endfor
 return,integ
end

;pro time_chihel, path, d3, hel

RR = read_major_radius(path)
path = './'
d3 = 0
hel = -10.
;#; understand which kind of energy is being computed'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
 mf_field='dat/imf_bprofiles.sav'
 dum = file_test(path+mf_field)
 print,'cacca',path+mf_field
 if (dum ne 1) then begin
  call_procedure, 'isavemodeprofiles',path,1,d3
 endif
 cmpt= test_up_to_date2(path,path+mf_field)
 print,'cmpt',cmpt
 if (cmpt gt 0) then begin
  call_procedure, 'isavemodeprofiles',path,1,d3
 endif
 
 restore, path+mf_field
 print, 'restored!'
 print, 'helicity=', hel

 nr1 = n_elements(pror1)
 nr2 = n_elements(pror2)
 ntt = n_elements(br[*,0,0,0])
 outfile = path+'/dat/rot_freq.sav'

;#; creo vettore con profilo temporale e radiale della rotazione del modo m=1
 adv = make_array(n_elements(br[*,0,0,0]),n_elements(br[0,*,0,0]),/double,value=0.)
 freq_rot = make_array(n_elements(br[*,0,0,0]),n_elements(br[0,*,0,0]),/double,value=0.)
 vel_rot = make_array(n_elements(br[*,0,0,0]),n_elements(br[0,*,0,0]),/double,value=0.)

 jmn = 1
 adv[*,*] = reform(br[*,*,jmn,1])
;#; prima di tutto raddrizzo la fase
 for i=0,nr1-1 do begin
  cnt = 0
  for q=0, ntt-2 do begin
   if ( (br[q+1,i,jmn,1] - br[q,i,jmn,1]) lt 0. ) then begin
    cnt = cnt + 1
;     print,i,q,cnt
    adv(q+1:n_elements(br[*,0,0,0])-1,i) = br(q+1:n_elements(br[*,0,0,0])-1,i,jmn,1) + 2. * !Dpi * cnt
;#; phase advances
   endif else begin
;#; phase goes back
   endelse
  endfor
  dum = deriv(tt,reform(adv[*,i]))
  freq_rot[*,i] = dum *  RR[0] / abs(hel)
  vel_rot[*,i] = freq_rot[*,i] * pror1[i]
 endfor
 
 nlev=20
 lev=findgen(nlev)/(nlev-1) * 0.05 +0.03
end
