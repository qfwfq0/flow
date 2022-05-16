pro antifft_guscio, path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp

;#; fa la antifft solo su un guscio a r=rgus

rgus = -1

;path = './'
RR = read_major_radius(path)
;d3 = 1
;hel = 10.
;str = 'j'
help,str
print,'field=',str
;#; understand which kind of quantity is being computed'
case str of
 'b': mssg='Reconstructing b field'
 'v': mssg='Reconstructing v field'
 'j': mssg='Reconstructing j field'
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


;restore,path+'dat/spectrum.sav'
;spectrum=read_spectrum(path)
;mm=spectrum.mm
;nanf=spectrum.nanf
;nz=spectrum.nz
;nzcon=spectrum.nzcon
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
 nzcon=nz
endelse
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
if (str eq 'j') then begin
 if (d3 eq 1) then begin
  if (n_elements(pjr[0,0,*]) lt total(nz[0:1])) then begin
   my = 1
  endif else begin
   my = 2
  endelse
 endif else begin
 endelse
endif
if (str eq 'b') then begin
 if (d3 eq 1) then begin
  if (n_elements(pdbr[0,0,*]) lt total(nz[0:1])) then begin
   my=1
  endif else begin
   my = 2
   print, 'my = ',my
  endelse
 endif else begin
 endelse
endif
if (str eq 'v') then begin
 if (d3 eq 1) then begin
  if (n_elements(pdvr[0,0,*]) lt total(nz[0:1])) then begin
   my=1
  endif else begin
   my = 2
  endelse
 endif else begin
 endelse
endif
my=2

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
  outfile = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'_guscio_m1n11-13.sav'
print,outfile
   print,'computing antifft of the '+str+' at time'+strtrim(string(pdt(itp(q)),format='(f11.4)'))+' time step'+strtrim(string(itp(q),format='(i0)'))

   if (str eq 'j') then begin
    antifft_r = make_array(nth,nzz,/double,value=0.)
    antifft_t = make_array(nth,nzz,/double,value=0.)
    antifft_z = make_array(nth,nzz,/double,value=0.)
   endif 
   if (str eq 'b') then begin
    antifft_r = make_array(nth,nzz,/double,value=0.)
    antifft_t = make_array(nth,nzz,/double,value=0.)
    antifft_z = make_array(nth,nzz,/double,value=0.)
   endif 
   if (str eq 'v') then begin
    antifft_r = make_array(nth,nzz,/double,value=0.)
   ; antifft_t = make_array(nitp,n_elements(pdvt[0,*,0]),nth,nzz,/double,value=0.)
   ; antifft_z = make_array(nitp,n_elements(pdvz[0,*,0]),nth,nzz,/double,value=0.)
   endif 

   print,'check',max(antifft_r) 
   for z = 0 , nzz - 1 do begin
    if (z mod 32 eq 0) then print,'z',z
    for t = 0 , nth - 1 do begin
;     rchi(*,t,z) = -rpdaz_0 + dieli / RR[0] * rpdat_0
      for im = 0, my do begin
       for in = nanf(im),nanf(im)+nzcon(im)-1 do begin
;#; Marco, per calcolare il caso m=2
        jmn = janz2(im,in,my,mm,nanf,nzcon)

;#; Marco, mode selection
;        if (q eq 0 and z eq 1 and t eq 1) then begin
;        print,'aaa',im,in,jmn
;        endif
;#; mode selection
;       if (im eq 0) then begin
;;         if (in ne 0) then goto, jump
;        endif
;       if (im eq 1) then begin
;        if ((in lt -25) or (in gt -1)) then goto, jump
;       endif

         if (im gt 2) then goto, jump
;#; solo modi m=0
;         if (im ge 1) then goto, jump
;#; solo modi m=2
;         if (im ne 2) then goto, jump

;#; solo parte m=1 e no assialsimmetric
;        if (im eq 0) then goto, jump
;        if (im ge 2) then goto, jump
;        if (im eq 1) then begin
;         if ((in lt -12) or (in gt -3)) then goto, jump
;        endif

;#; Marco, per tenere solo i due modi
;        if (im ne 1) then goto, jump
;        if (not ((in ne -11) xor (in ne -13))) then goto, jump

;#; Marco, per tenere solo i tre modi
;        if (im ne 1) then goto, jump
;        if ((in ne -7) xor (in ne -9) xor (in ne -11)) then goto, jump
;        if (not ((in ne -7) xor (in ne -10) xor (in ne -11) xor (in ne -13))) then goto, jump

        if (z eq 0 and t eq 0) then begin
         print,'modi considerati',im,in,jmn
        endif

;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro!!
        fcos =  cos(im * vecth(t) + in * vecz(z) / RR[0])
        fsin = - sin(im * vecth(t) + in * vecz(z) / RR[0])

         if (str eq 'j') then begin
          if (jmn eq 0) then begin
           antifft_r[t,z] = antifft_r[t,z] + real_part(pjr[itp[q],rgus,jmn])
           antifft_t[t,z] = antifft_t[t,z] + real_part(pjt[itp[q],rgus,jmn])
           antifft_z[t,z] = antifft_z[t,z] + real_part(pjz[itp[q],rgus,jmn])
          endif else begin
           antifft_r[t,z] = antifft_r[t,z] + real_part(reform(pjr[itp[q],rgus,jmn])) * fcos - imaginary(reform(pjr[itp[q],rgus,jmn])) * fsin
           antifft_t[t,z] = antifft_t[t,z] + real_part(reform(pjt[itp[q],rgus,jmn])) * fcos - imaginary(reform(pjt[itp[q],rgus,jmn])) * fsin
           antifft_z[t,z] = antifft_z[t,z] + real_part(reform(pjz[itp[q],rgus,jmn])) * fcos - imaginary(reform(pjz[itp[q],rgus,jmn])) * fsin
          endelse
         endif
         if (str eq 'b') then begin
          if (jmn eq 0) then begin
           antifft_r[t,z] = antifft_r[t,z] + real_part(pdbr[itp[q],rgus,jmn])
           antifft_t[t,z] = antifft_t[t,z] + real_part(pdbt[itp[q],rgus,jmn])
           antifft_z[t,z] = antifft_z[t,z] + real_part(pdbz[itp[q],rgus,jmn])
          endif else begin
           antifft_r[t,z] = antifft_r[t,z] + real_part(reform(pdbr[itp[q],rgus,jmn])) * fcos + imaginary(reform(pdbr[itp[q],rgus,jmn])) * fsin
           antifft_t[t,z] = antifft_t[t,z] + real_part(reform(pdbt[itp[q],rgus,jmn])) * fcos + imaginary(reform(pdbt[itp[q],rgus,jmn])) * fsin
           if (q eq 0 and z eq 1 and t eq 1) then begin
;            print,im,in,antifft_t[q,0:1,t,z], '   ', pdbt[itp[q],0:1,jmn]

           endif

           antifft_z[t,z] = antifft_z[t,z] + real_part(reform(pdbz[itp[q],rgus,jmn])) * fcos + imaginary(reform(pdbz[itp[q],rgus,jmn])) * fsin
          endelse
         endif
         if (str eq 'v') then begin
          if (jmn eq 0) then begin
           antifft_r[t,z] = antifft_r[t,z] + real_part(pdvr[itp[q],rgus,jmn])
;           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(pdvt[itp[q],*,jmn])
;           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(pdvz[itp[q],*,jmn])
          endif else begin
           antifft_r[t,z] = antifft_r[t,z] + real_part(reform(pdvr[itp[q],rgus,jmn])) * fcos - imaginary(reform(pdvr[itp[q],rgus,jmn])) * fsin
;           antifft_t[*,t,z] = antifft_t[*,t,z] + real_part(reform(pdvt[itp[q],*,jmn])) * fcos - imaginary(reform(pdvt[itp[q],*,jmn])) * fsin
;           antifft_z[*,t,z] = antifft_z[*,t,z] + real_part(reform(pdvz[itp[q],*,jmn])) * fcos - imaginary(reform(pdvz[itp[q],*,jmn])) * fsin
          endelse
         endif

        jump: 
       endfor
      endfor
    endfor
   endfor

  ptime=pdt(itp)
  print,outfile
  if (str eq 'b') then begin
   br = antifft_r
   bt = antifft_t
   bz = antifft_z
   save,filename=outfile,br,bt,bz,itp,pror1,pror2,vecth,vecz,ptime
;   save,filename=outfile,br,itp,pror1,pror2,vecth,vecz,ptime
  endif
  if (str eq 'v') then begin
   vr = antifft_r
   vt = antifft_t
   vz = antifft_z
   save,filename=outfile,vr,vt,vz,itp,pror1,pror2,vecth,vecz,ptime
  endif
  if (str eq 'j') then begin
   jr = antifft_r
   jt = antifft_t
   jz = antifft_z
   save,filename=outfile,jr,jt,jz,itp,pror1,pror2,vecth,vecz,ptime
;   save,filename=outfile,jr,itp,pror1,pror2,vecth,vecz,ptime
  endif
  endfor ; end of cycle on time

end
