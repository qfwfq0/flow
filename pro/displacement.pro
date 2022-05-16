pro displacement, path, d3,nth,nzz,itp0, itp1, deltaitp

;#; Marco, computes the "radial displacement" associated to one MHD mode or to a set of MHD modes, based on the formula csi_r^(mn) = (r*br^(mn)/(i*Bt00*(m+nq)))
; path = './'
 RR = read_major_radius(path)
print,'mesh',nth,nzz
; d3 = 1
; hel = 10.
; str = 'b'
; itp0 = 1900
; itp1 = 1900
; deltaitp=1
;#; understand which kind of quantity is being computed'
 mf_field='dat/ibprofiles.sav'
 dum = file_test(path+mf_field)
 if (dum ne 1) then begin
  print,'need to create the fields file ibprofiles.sav',d3
  call_procedure,'isavemodeprofiles',path,'0',d3
;       print,'implement this part'
;       stop
;       call_procedure,'imergebmf',path
 endif

 if (pdbr eq !null) then begin
  restore, path+mf_field
 endif

result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
;spectrum=read_spectrum(path)
;mm=spectrum.mm
;nanf=spectrum.nanf
;nz=spectrum.nz
;nzcon=spectrum.nzcon
;itp0=100
;itp1=105
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


;#; check on saved modes
 if (d3 eq 1) then begin
   my=1
   im = 1
   nmin = -15
   nmax = -6
   jmin =  67
   jmax = 75
   j0min =  1
   j0max = 11
   nmodes = jmax - jmin
   nmodes0 = j0max - j0min
;   jmin = janz2(im,nmin,my,mm,nanf,nzcon)
;   jmax = janz2(im,nmax,my,mm,nanf,nzcon)
 endif else begin
   my = 1
   jmin = 1
   jmax = 9
 endelse

 unit_i = dcomplex(0.d,1.d)
 nr1 = n_elements(pror1)
 nr2 = n_elements(pror2)
; nth = 64
; nzz =  128
 vecth = dindgen(nth)/(nth-1) * 2. * !Dpi
 if (nzz eq 1) then begin
  vecz = 0.
 endif else begin
  vecz = dindgen(nzz)/(nzz-1) * 2. * !Dpi * RR[0]
 endelse

 qqq = make_array(n_elements(pdbt[0,*,0]),/double,value=0.d)
 for q = 0 , nsaved-1 do begin ;#; beginning of the cycle in time
 ;#; compute axisymmetric safety factor
  qqq = (pror2 * pdbz[itp(q),*,0]) / (RR(0) * pdbt[itp(q),*,0]) 
  print,'time=',q,itp(q)
  spawn, 'mkdir -p '+path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))
  outfile = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/displacement_m1.sav'
  print,outfile
  displm1 = make_array(nmodes+1,n_elements(pdbt[0,*,0]),nth,nzz,/dcomplex,value=0.d)
  tot_displm1 = make_array(n_elements(pdbt[0,*,0]),nth,nzz,/dcomplex,value=0.d)
  displm0 = make_array(nmodes0+1,n_elements(pdbt[0,*,0]),nzz,/dcomplex,value=0.d)
  tot_displm0 = make_array(n_elements(pdbt[0,*,0]),nzz,/dcomplex,value=0.d)
;  displm1 = make_array(nmodes+1,n_elements(pdbt[0,*,0]),nth,nzz,/double,value=0.d)

;#; section about m=1 modes
  for p = 0,nmodes do begin
;#; interpolate pdbr on mesh x2
   mnmn = mnum(p+jmin,nz,nanf,d3) 
   if (q eq 0) then begin
   print,p,p+jmin,mnmn
   endif
   dbr = interpol(reform(pdbr[itp(q),*,p+jmin]), pror1, pror2 )
   if (mnmn[0] ne 1) then begin
    dbr[0] = - dbr[1]
   endif else begin
    dbr[0] = dbr[1]
   endelse
;#; start reconstruction cycle
   for z = 0 , nzz - 1 do begin
    for t = 0 , nth - 1 do begin

     facty = (pror2 * dbr) / (unit_i * (mnmn(0) + mnmn(1)*qqq) * pdbt[itp(q),*,0]) * exp(unit_i*(mnmn(0)*vecth(t) + mnmn(1)/RR(0)*vecz(z)))
     if (p eq 7) then begin 
;     plot,pror2,facty,title='m='+strtrim(string(mnmn(0),'(i0)'))+' n='+strtrim(string(mnmn(1),'(i0)'))
;pait, 0.2
     endif
     displm1(p,*,t,z) = facty
;     displm1(p,*,t,z) = abs(facty)

    endfor ;#; end of antifft cycle in t
   endfor ;#; end of antifft cycle in z
  endfor ;#; end of cycle in modes

;#; section about m=0 modes, no theta dependence
  for p = 0,nmodes0 do begin
;#; interpolate pdbr on mesh x2
   mnmn = mnum(p+j0min,nz,nanf,d3) 
;   if (q eq 0) then begin
;   endif
;   dbr = interpol(reform(pdbr[itp(q),*,p+j0min]), pror1, pror2 )
;   if (mnmn[0] ne 1) then begin
;    dbr[0] = - dbr[1]
;   endif else begin
;    dbr[0] = dbr[1]
;   endelse
;#; start reconstruction cycle
   for z = 0 , nzz - 1 do begin

     facty = (pror2 * dbr) / (unit_i * (mnmn(0) + mnmn(1)*qqq) * pdbt[itp(q),*,0]) * exp(unit_i*(mnmn(1)/RR(0)*vecz(z)))
;     if (p eq 7) then begin 
;     plot,pror2,facty,title='m='+strtrim(string(mnmn(0),'(i0)'))+' n='+strtrim(string(mnmn(1),'(i0)'))
;pait, 0.2
;     endif
     displm0(p,*,z) = facty
;     displm1(p,*,t,z) = abs(facty)

   endfor ;#; end of antifft cycle in z
  endfor ;#; end of cycle in modes

;#; linear sum of displacements associated with m=0 modes
 tot_displm0 = total(displm0,1)
;#; linear sum of displacements associated with m=1 modes
 tot_displm1 = total(displm1,1)
;#; save part
 save,filename=outfile,itp,jmin,jmax,pror1,pror2,vecth,vecz,displm1,tot_displm1,displm0,tot_displm0
 endfor ;#; end of the cycle in time
end

