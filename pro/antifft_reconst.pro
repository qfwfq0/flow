pro antifft_reconst, path, str, d3, hel

RR = read_major_radius(path)
;path = './'
;d3 = 0
;hel = 10.
;str = 'b'
;#; understand which kind of energy is being computed'
case str of
 'b': mssg='Reconstructing b field'
 'v': mssg='Reconstructing v field'
 'j': mssg='Reconstructing j field'
endcase
case str of
 'b': begin
      mf_field='dat/imf_bprofiles.sav'
      dum = file_test(path+mf_field)
      if (dum ne 1) then begin
       call_procedure,'imergebmf',path
      endif
      end
 'v': begin
      mf_field='dat/imf_vprofiles.sav'
      dum = file_test(path+mf_field)
      if (dum ne 1) then begin
       call_procedure,'imergevmf',path
      endif
      end
 'j': begin
      mf_field='dat/imf_jprofiles.sav'
      dum = file_test(path+mf_field)
      if (dum ne 1) then begin
       mf = 1
       call_procedure,'isavejprofiles',path,mf,d3
      endif
      end
endcase

restore, path+mf_field
print, 'restored!'

nr1 = n_elements(pror1)
nr2 = n_elements(pror2)
nth = 32
nz =  8
vectth = dindgen(nth)/(nth-1) * 2. * !Dpi
vecz = dindgen(nz)/(nz-1) * 2. * !Dpi * RR[0]

outfile = path+'/dat/antifft_'+str+'.sav'
if (d3 eq 0) then begin
 dum = file_test(outfile)
 tme_sav = file_info(outfile)
 time_sav = tme_sav.mtime 
 tme_profsav = file_info(path+mf_field)
 time_profsav = tme_profsav.mtime 
 if (dum ne 1 or time_sav lt time_profsav) then begin
 print,'2d reconstruction with hel= ',hel
 maxm = 4
 bbr = make_array(n_elements(br[*,0,0,0]),nr1,nth,nz)
;#; cycle on time
 for q = 0, n_elements(br[*,0,0,0]) - 1 do begin
  if (q mod 10 eq 0) then print,'computing itp= ',q
;#; cycle on mode number
  for mm = 0, maxm do begin
;#; cycle on nth
   for th = 0, nth-1 do begin
;#; cycle on nz
    for z = 0, nz-1 do begin
     bbr[q,*,th,z] = bbr[q,*,th,z] + br[q,*,mm,0] * real_part( exp( complex(0.,mm*vectth(th) - hel/RR[0] * vecz(z) + br[q,*,mm,1])))
    endfor
   endfor
  endfor
 endfor

;#; save
 save,filename=outfile,bbr,pror1,pror2,vectth,vecz
 endif else begin
;#; restore already computed
 restore,filename=outfile
 endelse
endif else begin
 print,'3d not implemented yet'
endelse



loadct,'0',file='~/idl/colors.tbl'
col=findgen(n_elements(br[*,0,0,0]))/(n_elements(br[*,0,0,0]) - 1)*253.
;window,1
; for q = 0, n_elements(br[*,0,0,0]) - 1 do begin
;    contour,reform(bbr[q,*,*,0]),pror1,vectth,nlev=20,color=col(q)
;   wait,0.02
; endfor
window,1
 for q = 0, n_elements(br[*,0,0,0]) - 1 do begin
    contour,reform(bbr[q,*,*,0]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=20,color=col(q),/fill,title=string(q)
   wait,0.12
 endfor
;window,2
; for q = 0, n_elements(br[*,0,0,0]) - 1 do begin
;    contour,reform(bbr[q,*,0,*]),pror1,vecz,nlev=20,color=col(q),/fill,title=string(q)
;   wait,0.12
; endfor

end

