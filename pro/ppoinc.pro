;pro ppoincs, path, itp, d3, hel

RR = read_major_radius(path)
path = './'
d3 = 1
hel = -8.
itp = '30'


nr = 0
nmodes = 0
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
binfile = path+'dat/itp/'+itp+'/bcyl.bin'
openr,lun,binfile,/get_lun,/f77_unformatted
 readu, lun, nr
 readu, lun, nmodes
 print,nr,nmodes
 br = make_array(nr,nmodes,/dcomplex)
 bt = make_array(nr+1,nmodes,/dcomplex)
 bz = make_array(nr+1,nmodes,/dcomplex)
 readu, lun, br
 readu, lun, bt
 readu, lun, bz
 
 free_lun,lun
print, 'read!'

;#; Marco, general data for NEMATO input
bcond = make_array(6,/integer)
bcond(0) = 7
bcond(1) = 0
bcond(2) = 1
bcond(3) = 1
bcond(4) = 1
bcond(5) = 1
tor = 1

nth = 32
nzz =  8
pror1 = dindgen(nr)/(nr-1) 
pror2 = dindgen(nr+1)/(nr-1) - 0.5d/(nr-1)
vectth = dindgen(nth)/(nth-1) * 2. * !Dpi
vecz = dindgen(nzz)/(nzz-1) * 2. * !Dpi * RR[0]

;#; Marco, Poincaré file name

dir = path+'/dat/itp/'+itp+'/nemato/'
call_procedure,'check_dir',dir

outfile = path+'/dat/itp/'+itp+'/nemato/poinc'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'.bin'
 openw,lun,outfile,/get_lun,/f77_unformatted

if (d3 eq 0) then begin
 dum = file_test(outfile)
 print,'2d reconstruction with hel= ',hel
; bbr = make_array(n_elements(br[*,0,0,0]),nr1,nth,nz)
;;#; cycle on time
; for q = 0, n_elements(br[0,*]) - 1 do begin
;  aaa=mnum(q,nz,nanf,d3)
;function mnum,j,nz,nanf,d3
;;#; cycle on mode number
;  for mm = 0, maxm do begin
;;#; cycle on nth
;   for th = 0, nth-1 do begin
;;#; cycle on nz
;    for z = 0, nz-1 do begin
;     bbr[q,*,th,z] = bbr[q,*,th,z] + br[q,*,mm,0] * real_part( exp( complex(0.,mm*vectth(th) - hel/RR[0] * vecz(z) + br[q,*,mm,1])))
;    endfor
;   endfor
;  endfor

;#; save
; save,filename=outfile,bbr,pror1,pror2,vectth,vecz

;#; 3D case
endif else begin
 
 bbr = make_array(nr+1,nth,nzz,/double,value=0.)
 bbt = make_array(nr+1,nth,nzz,/double,value=0.)
 bbz = make_array(nr+1,nth,nzz,/double,value=0.)

;#; first: interpolation of br on the x2 mesh
  ddbr =  bt * dcomplex(0.d,0.d)
  help,bt
  print,nanf
   for q=0, n_elements(bt[*,0])-1 do begin
   for j=0, n_elements(bt[0,*])-1 do begin
         dbr = interpol(reform(br[*,j]), pror1, pror2 )
         mn = mnum(j,nz,nanf,d3)
         if (mn[0] eq 0) then begin
          dbr[0] = -dbr[1]
         endif else begin
          if (mn[0] eq 1) then begin
           dbr[0] = dbr[1]
          endif else begin
           dbr[0] = -dbr[1]
          endelse
         endelse
     ddbr[*,j] = dbr
   endfor
   endfor

 for q = 0, n_elements(br[0,*]) - 1 do begin
  mn=mnum(q,nz,nanf,d3)
  print,'dealing with mode',mn
;#; cycle on nth
   for th = 0, nth-1 do begin
;#; cycle on nz
    for z = 0, nzz-1 do begin
;#; cycle on nz
     bbr[*,th,z] = bbr[*,th,z] + real_part(ddbr[*,q]) * 2. * cos(mn[0] * vectth(th) - mn[1] / RR[0] * vecz(z)) + imaginary(ddbr[*,q]) * 2. * sin(mn[0] * vectth(th) - mn[1] / RR[0] * vecz(z))
     bbt[*,th,z] = bbt[*,th,z] + real_part(bt[*,q]) * 2. * cos(mn[0] * vectth(th) - mn[1] / RR[0] * vecz(z)) + imaginary(bt[*,q]) * 2. * sin(mn[0] * vectth(th) - mn[1] / RR[0] * vecz(z))
     bbz[*,th,z] = bbz[*,th,z] + real_part(bz[*,q]) * 2. * cos(mn[0] * vectth(th) - mn[1] / RR[0] * vecz(z)) + imaginary(bz[*,q]) * 2. * sin(mn[0] * vectth(th) - mn[1] / RR[0] * vecz(z))
    endfor
;#; cycle on nth
   endfor

 endfor
      nrr = fix(nr+1)

      writeu, lun, nrr,nth,nzz,bcond,tor
      writeu, lun, pror2
      writeu, lun, vectth
      writeu, lun, vecz
      writeu, lun, bbr
      writeu, lun, bbt
      writeu, lun, bbz

endelse

free_lun,lun
end
