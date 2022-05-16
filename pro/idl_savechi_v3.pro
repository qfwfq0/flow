function integralc,vect,xx
;#; general purpose function to calculate the integral of a function (vect) along the mesh x2 (vect)
;#; x2=findgen(lx+2)/lx-1./(2.*lx)
;#; the result will be stored in the mesh x1(i) = findgen(lx+1)/lx

 integ = make_array(n_elements(xx)-1,/dcomplex,value=0.)
 yy = findgen(n_elements(xx)-1) / (n_elements(xx)-2)
 dx = 1.d / (n_elements(xx)-2)
 
; for q=0, n_elements(vect(*))-1 do begin
  integ(0) = 0.d
  for i=1,n_elements(yy)-1 do begin
;#; marco, 5 marzo 2019
;#; attenzione, la mesh IDL non parte da -1 come quella fortran quindi bisogna usare vect(i) invece che vect(i-1)
   integ(i) = integ(i-1) + vect(i) * dx
  endfor
; endfor
 return,integ
end

pro idl_savechi_v3,path,d3,tok,ieli,nth,nzz,itp0,itp1,deltaitp
;#; this routine computes the helical flux function.
;#; it can be computed on a 3D grids or on a 2D grid (both at th=cost,phi=cost)
;#; Marco Veranda, febbraio-marzo 2019

 RR = read_major_radius(path)

 result = file_test(path+'dat/spectrum.sav')
 if (result ne 1) then begin
  call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
  if (m2d eq 0) then begin
   m2d = 1./ieli
  endif
 endif else begin
  restore,path+'dat/spectrum.sav'
  if (nzcon eq !null) then begin
  nzcon=nz
  nzcon[0]=nz[0]+1
  endif
  if (m2d eq !null) then begin
   if (d3 eq 0) then begin
    m2d = 1
   endif else begin
    m2d = 0
   endelse
  endif
 endelse
 print,'ieli=',ieli,' , itp0=',itp0,'  m2d=',m2d,'   R0=',RR[0]
 if (nzcon[0] = nz[0]) then begin 
  nzcon[0] = nzcon[0] + 1
 endif
 print,'AAA',nzcon

;#; inputfile for magnetic field
 input_b_file = path+'dat/ibprofiles.sav'
 exx = file_test(input_b_file) 
 if (exx ne 1) then begin
  call_procedure, 'isavemodeprofiles',path,0,d3
  restore, input_b_file
 endif else begin
 ;#; check whether ibprofiles is up to date
  test = test_up_to_date_files(input_b_file,path+'dat/imf_bprofiles.sav')
  if (test eq 1) then begin
   call_procedure, 'isavemodeprofiles',path,0,d3
   restore, input_b_file
  endif else begin
   restore, input_b_file
  endelse 
 endelse

;#; create necessary directories
; if (it eq -1) then begin
;  it = n_elements(pdbr[*,0,0])-1
; endif

;#; Marco, deal with time
 if (itp0 eq !null) then begin
  itp0 = 0
 endif
 if (itp1 eq !null) then begin
  itp1 = n_elements(pdt) - 1
 endif
 
 if (itp0 eq !null) then begin
  print,'define itp0!'
 endif
 itp0 = fix(itp0,type=3)
 itp1 = fix(itp1,type=3)
 nitp = itp1 - itp0 + 1
 nsaved = ceil(float(nitp)/deltaitp)
 itp= indgen(nsaved)*deltaitp + itp0 


 nr = n_elements(pdbr[0,*,0])
 pror1 = dindgen(nr) / (nr - 1)
 th = dindgen(nth) / (nth - 1) * 2. * !Dpi
 zz = dindgen(nzz) / (nzz - 1) * 2. * !Dpi * RR[0]
 rchi_tot = make_array(nsaved,nr,nth,nzz,/double,value=0.)

 if (tok eq 1) then begin
  maxm = 9
 endif else begin
  maxm = max(mm)
  if (maxm gt 4) then begin
   maxm = 4
  endif
 endelse 
 ii = dcomplex(0.d,1.d)
 dieli = double(ieli)

;#; ciclo sui tempi
 for q=0,nsaved-1 do begin
 rchi = make_array(nr,nth,nzz,/double,value=0.)
 call_procedure, 'check_dir',path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))
;#; intanto integro la parte assialsimmetrica del campo per creare Az00, e rAtheta00
 az00 = make_array(n_elements(pdbr[0,*,0]),/dcomplex,value=0.)
 rat00 = make_array(n_elements(pdbr[0,*,0]),/dcomplex,value=0.)
 rat00 = integralc( reform( pdbz[itp(q),*,0] ) * pror2,pror2 )
 az00 = integralc( - reform( pdbt[itp(q),*,0] ) ,pror2 )
; help, az00
; help, rat00

; if (tok eq 1) then begin
 for z = 0 , nzz - 1 do begin
  if (z mod 32 eq 0) then print,'z',z,' dieli',dieli
  for t = 0 , nth - 1 do begin
   rchi(*,t,z) = az00 - dieli/RR[0] * rat00

   for im = 1, maxm do begin
;#; spectrum selection part
    if (d3 eq 0) then begin
     mmf = im * m2d
     jmn = im
    endif else begin
     mmf = im
     jmn = janz2(im,dieli*im,my,mm,nanf,nzcon)
    endelse
;     print,'m=',mmf,'  n=',dieli*mmf,' jmn=',jmn,'  itp(q)=',itp(q)
;     print, 'casino',jmn,floor(jmn),abs(floor(jmn)-jmn)
    if (abs(floor(jmn)-jmn) gt 0.01) then begin
;;     print, 'casinovero',jmn,floor(jmn),abs(floor(jmn)-jmn)
     continue
    endif
    if (z eq 0 and t eq 0) then begin
     print,'m=',mmf,'  n=',dieli*mmf,' jmn=',jmn,'  itp(q)=',itp(q)
    endif
;;#; in ibprofiles.sav c'è già il fattore 2, quindi qui non lo metto
;     rchi(*,t,z) = rchi(*,t,z) + pror1 / mmf * (real_part(pdbr[itp(q),*,jmn]) * cos(mmf * th(t) + mmf * dieli * zz(z) / RR[0] - !pi/2.) - imaginary(pdbr[itp(q),*,jmn]) * sin(mmf * th(t) + mmf * dieli * zz(z) / RR[0] - !pi/2.))
     rchi(*,t,z) = rchi(*,t,z) + pror1 / mmf * (real_part(pdbr[itp(q),*,jmn]) * sin(mmf * th(t) + mmf * dieli * zz(z) / RR[0]) + imaginary(pdbr[itp(q),*,jmn]) * cos(mmf * th(t) + mmf * dieli * zz(z) / RR[0] ))
   endfor

  endfor
 endfor

 rchi_tot(q,*,*,*) = rchi
; endif else begin ;#; begin RFP part
; endelse
 inteli=fix(dieli)
 save_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chi'+strtrim(string(inteli,format='(i0)'))+'.sav'
 print,'saving in.... ',save_file
 save,filename=save_file,rchi,pror1,th,zz

 call_procedure,'plot_chihel', path, d3, tok, dieli,itp(q)
 endfor ;#; end of cycle on times

 save_file_tot = path+'dat/itp/chitot_ieli'+strtrim(string(inteli,format='(i0)'))+'_'+strtrim(string(itp0,format='(i0)'))+'-'+strtrim(string(itp1,format='(i0)'))+'-'+strtrim(string(deltaitp,format='(i0)'))+'.sav'
 save,filename=save_file_tot,rchi_tot,pror1,th,zz,itp,nsaved,dieli

  call_procedure, 'show_time_chihel',path,itp0,itp1,inteli,deltaitp,itp
end
