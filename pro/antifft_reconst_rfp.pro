pro antifft_reconst_rfp, path, str, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax

;#; CHANGELOG
;#; Marco, 10 dic 2019: avevo sbagliato un segno nella antitrasformata di Fourier

;path = './'
print,'aaaa',nth
RR = read_major_radius(path)
;d3 = 1
;hel = 10.
;str = 'j'
help,str
print,'field=',str
mymin=fix(mymin)
mymax=fix(mymax)
;#; understand which kind of quantity is being computed'
case str of
 'b': mssg='Reconstructing b field'
 'v': mssg='Reconstructing v field'
 'j': mssg='Reconstructing j field'
 'om': mssg='Reconstructing vorticity field'
endcase
case str of
 'b': begin
      mf_field='dat/ibprofiles.sav'
      dum = file_test(path+mf_field)
      if (dum ne 1) then begin
       print,'need to create the fields file ibprofiles.sav',d3
       call_procedure,'isavemodeprofiles',path,'0',d3
;       print,'implement this part'
;       stop
;       call_procedure,'imergebmf',path
      endif
      end
 'v': begin
      mf_field='dat/ivprofiles.sav'
      dum = file_test(path+mf_field)
      if (dum ne 1) then begin
       print,'need to create the fields file ibprofiles.sav',d3
       call_procedure,'isavemodeprofiles',path,'0',d3
;       print,'implement this part'
;       stop
;       call_procedure,'imergevmf',path
      endif
      end
 'j': begin
      mf_field='dat/ijprofiles.sav'
      dum = file_test(path+mf_field)
      if (dum ne 1) then begin
       mf = 0
       call_procedure,'isavejprofiles',path,mf,d3
      endif
      end
 'om': begin
      mf_field='dat/iomprofiles.sav'
      dum = file_test(path+mf_field)
      if (dum ne 1) then begin
       mf = 0
       call_procedure,'isaveomprofiles',path,mf,d3
      endif
      end
endcase

if (str eq 'b' and pdbr eq !null) then begin
 restore, path+mf_field
endif
if (str eq 'v' and pdbr eq !null) then begin
 restore, path+mf_field
endif
if (str eq 'j' and pjr eq !null) then begin
 restore, path+mf_field
endif
if (str eq 'om' and pomr eq !null) then begin
 restore, path+mf_field
endif


;restore,path+'dat/spectrum.sav'
;spectrum=read_spectrum(path)
;mm=spectrum.mm
;nanf=spectrum.nanf
;nz=spectrum.nz
;nzcon=spectrum.nzcon


result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
 print,'nzcon',nzcon
endif else begin
 restore,path+'dat/spectrum.sav'
 if (nzcon eq !null) then begin
  nzcon = nz
  if ((nzcon[0]-nz[0]) gt 1 ) then stop
  print,nzcon
 endif
endelse

if (nzcon[0] = nz[0]) then begin
 nzcon[0] = nzcon[0] + 1
endif
  print,'antifft_reconst_rfp, nzcon: ',nzcon
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
if (itp1 eq -1) then begin
itp1 = n_elements(pdt) - 1
endif
if (itp0 eq -1) then begin
itp0 = n_elements(pdt) - 1
endif
itp0 = fix(itp0,type=3)
itp1 = fix(itp1,type=3)
nitp = itp1 - itp0 + 1
nsaved = ceil(float(nitp)/deltaitp)
itp= indgen(nsaved)*deltaitp + itp0 

;#; check on saved modes
if (str eq 'j') then begin
 if (d3 eq 1) then begin
  if (n_elements(pjr[0,0,*]) lt total(nz[0:1])) then begin
   my = mymax
  endif else begin
   my = mymax
  endelse
 endif else begin
   my = mymax
 endelse
endif
if (str eq 'b') then begin
 if (d3 eq 1) then begin
  if (n_elements(pdbr[0,0,*]) lt total(nz[0:1])) then begin
   my=mymax
  endif else begin
   my = mymax
   print, 'my = ',my
  endelse
 endif else begin
   my = mymax
 endelse
endif
if (str eq 'v') then begin
 if (d3 eq 1) then begin
  if (n_elements(pdvr[0,0,*]) lt total(nz[0:1])) then begin
   my = mymax
  endif else begin
   my = mymax
  endelse
 endif else begin
   my = mymax
 endelse
endif
if (str eq 'om') then begin
 if (d3 eq 1) then begin
  if (n_elements(pomr[0,0,*]) lt total(nz[0:1])) then begin
   my = mymax
  endif else begin
   my = mymax
  endelse
 endif else begin
   my = mymax
 endelse
endif

nr1 = n_elements(pror1)
nr2 = n_elements(pror2)
;nth = 128
;nzz =  1
vecth = dindgen(nth)/(nth-1) * 2. * !Dpi
if (nzz eq 1) then begin
 vecz = 0.
endif else begin
 vecz = dindgen(nzz)/(nzz-1) * 2. * !Dpi * RR[0]
endelse

;#; check directory existence

;#; begin 3d part
  for q = 0 , nsaved-1 do begin
  print,'time=',q,itp(q)
  spawn, 'mkdir -p '+path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))
  outfile = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(102,format='(i0)'))+'x'+strtrim(string(nth,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
print,outfile
   print,'computing antifft of the '+str+' at time'+strtrim(string(pdt(itp(q)),format='(f11.4)'))+' time step'+strtrim(string(itp(q),format='(i0)'))

   if (str eq 'j') then begin
    antifft_r = make_array(n_elements(pjr[0,*,0]),nth,nzz,/double,value=0.)
    antifft_t = make_array(n_elements(pjt[0,*,0]),nth,nzz,/double,value=0.)
    antifft_z = make_array(n_elements(pjz[0,*,0]),nth,nzz,/double,value=0.)
   endif 
   if (str eq 'om') then begin
    antifft_r = make_array(n_elements(pomr[0,*,0]),nth,nzz,/double,value=0.)
    antifft_t = make_array(n_elements(pomt[0,*,0]),nth,nzz,/double,value=0.)
    antifft_z = make_array(n_elements(pomz[0,*,0]),nth,nzz,/double,value=0.)
   endif 
   if (str eq 'b') then begin
    antifft_r = make_array(n_elements(pdbr[0,*,0]),nth,nzz,/double,value=0.)
    antifft_t = make_array(n_elements(pdbt[0,*,0]),nth,nzz,/double,value=0.)
    antifft_z = make_array(n_elements(pdbz[0,*,0]),nth,nzz,/double,value=0.)
   endif 
   if (str eq 'v') then begin
    antifft_r = make_array(n_elements(pdvr[0,*,0]),nth,nzz,/double,value=0.)
    antifft_t = make_array(n_elements(pdvt[0,*,0]),nth,nzz,/double,value=0.)
    antifft_z = make_array(n_elements(pdvz[0,*,0]),nth,nzz,/double,value=0.)
   endif 

   print,'check',max(antifft_r) 
   for z = 0 , nzz - 1 do begin
    if (z mod 32 eq 0) then print,'z',z
    for t = 0 , nth - 1 do begin
;     rchi(*,t,z) = -rpdaz_0 + dieli / RR[0] * rpdat_0
      for im = mymin, mymax do begin
       for in = nanf(im),nanf(im)+nzcon(im)-1 do begin
;        print,im,my,nanf(im),nanf(im)+nzcon(im)-1
        jmn = janz2(im,in,my,mm,nanf,nzcon)

;        if (q eq 0 and z eq 1 and t eq 1) then begin
;        print,'aaa',im,in,jmn
;        endif
;#; mode selection
;       if (im eq 0) then begin
;         if (in ne 0) then goto, jump
;        endif
;        if (im eq 1) then begin
;         if ((in lt -25) or (in gt -3)) then goto, jump
;        endif

;#; solo modi m=0
;         if (im ge 1) then goto, jump

;#; solo parte m=1 e no assialsimmetric
;        if (im eq 0) then goto, jump
;        if (im eq 1) then begin
;         if ((in lt -12) or (in gt -3)) then goto, jump
;        endif

        if (z eq 1 and t eq 1) then begin
         print,'modi considerati',im,in,jmn
        endif

;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro!!
        fcos =  cos(im * vecth(t) + in * vecz(z) / RR[0])
        fsin = - sin(im * vecth(t) + in * vecz(z) / RR[0])

         if (str eq 'j') then begin
          if (jmn eq 0) then begin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(pjr[itp[q],*,jmn])
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(pjt[itp[q],*,jmn])
           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(pjz[itp[q],*,jmn])
          endif else begin
           ;#; Marco, un segno meno è già dentro fsin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(reform(pjr[itp[q],*,jmn])) * fcos + imaginary(reform(pjr[itp[q],*,jmn])) * fsin
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(reform(pjt[itp[q],*,jmn])) * fcos + imaginary(reform(pjt[itp[q],*,jmn])) * fsin
           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(reform(pjz[itp[q],*,jmn])) * fcos + imaginary(reform(pjz[itp[q],*,jmn])) * fsin
          endelse
         endif
         if (str eq 'om') then begin
          if (jmn eq 0) then begin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(pomr[itp[q],*,jmn])
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(pomt[itp[q],*,jmn])
           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(pomz[itp[q],*,jmn])
          endif else begin
           ;#; Marco, un segno meno è già dentro fsin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(reform(pomr[itp[q],*,jmn])) * fcos + imaginary(reform(pomr[itp[q],*,jmn])) * fsin
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(reform(pomt[itp[q],*,jmn])) * fcos + imaginary(reform(pomt[itp[q],*,jmn])) * fsin
           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(reform(pomz[itp[q],*,jmn])) * fcos + imaginary(reform(pomz[itp[q],*,jmn])) * fsin
          endelse
         endif
         if (str eq 'b') then begin
          if (jmn eq 0) then begin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(pdbr[itp[q],*,jmn])
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(pdbt[itp[q],*,jmn])
           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(pdbz[itp[q],*,jmn])
          endif else begin
           ;#; Marco, un segno meno è già dentro fsin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(reform(pdbr[itp[q],*,jmn])) * fcos + imaginary(reform(pdbr[itp[q],*,jmn])) * fsin
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(reform(pdbt[itp[q],*,jmn])) * fcos + imaginary(reform(pdbt[itp[q],*,jmn])) * fsin
           if (q eq 0 and z eq 1 and t eq 1) then begin
;            print,im,in,antifft_t[q,0:1,t,z], '   ', pdbt[itp[q],0:1,jmn]

           endif

           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(reform(pdbz[itp[q],*,jmn])) * fcos + imaginary(reform(pdbz[itp[q],*,jmn])) * fsin
          endelse
         endif
         if (str eq 'v') then begin
          if (jmn eq 0) then begin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(pdvr[itp[q],*,jmn])
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(pdvt[itp[q],*,jmn])
           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(pdvz[itp[q],*,jmn])
          endif else begin
           ;#; Marco, un segno meno è già dentro fsin
           antifft_r[*,t,z] = antifft_r[*,t,z] + real_part(reform(pdvr[itp[q],*,jmn])) * fcos + imaginary(reform(pdvr[itp[q],*,jmn])) * fsin
           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(reform(pdvt[itp[q],*,jmn])) * fcos + imaginary(reform(pdvt[itp[q],*,jmn])) * fsin
           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(reform(pdvz[itp[q],*,jmn])) * fcos + imaginary(reform(pdvz[itp[q],*,jmn])) * fsin
          endelse
         endif

        jump: 
       endfor
      endfor
    endfor
   endfor

  ptime=pdt(itp)
  print,'inside antifft_reconst_rfp',outfile
  if (str eq 'b') then begin
   br = reform(antifft_r)
   bt = reform(antifft_t)
   bz = reform(antifft_z)
   save,filename=outfile,br,bt,bz,itp,pror1,pror2,vecth,vecz,ptime
;   save,filename=outfile,br,itp,pror1,pror2,vecth,vecz,ptime
  endif
  if (str eq 'v') then begin
   vr = reform(antifft_r)
   vt = reform(antifft_t)
   vz = reform(antifft_z)
   save,filename=outfile,vr,vt,vz,itp,pror1,pror2,vecth,vecz,ptime
  endif
  if (str eq 'j') then begin
   jr = reform(antifft_r)
   jt = reform(antifft_t)
   jz = reform(antifft_z)
   save,filename=outfile,jr,jt,jz,itp,pror1,pror2,vecth,vecz,ptime
;   save,filename=outfile,jr,itp,pror1,pror2,vecth,vecz,ptime
  endif
  if (str eq 'om') then begin
   omr = reform(antifft_r)
   omt = reform(antifft_t)
   omz = reform(antifft_z)
   save,filename=outfile,omr,omt,omz,itp,pror1,pror2,vecth,vecz,ptime
;   save,filename=outfile,jr,itp,pror1,pror2,vecth,vecz,ptime
  endif
  endfor ; end of cycle on time

end
