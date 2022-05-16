pro i_el, path,d3,mf,equ
;#; program computing the electric field in cylindrical geometry

; path='./'
; path='/disk26a/ricercatori/ft/specyl/veranda/flow/archive/rfp/3d/theta16/s100m10/base/'
; d3=1
; equ = 1
; itp=200

;#; test wheter file already exists!
 savfile = path+'dat/imf_bprofiles.sav'
 if (equ eq 1) then begin
  eq_file = path+'dat/el_eq.sav'
  uptodate =test_up_to_date_files(eq_file,savfile)
  if (uptodate eq 0) then begin
   goto, endprogram
  endif
 endif else begin
  el_file = path+'dat/el.sav'
  uptodate =test_up_to_date_files(el_file,savfile)
  if (uptodate eq 0) then begin
   goto, endprogram
  endif
 endelse

; restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse

 nsectionsdat = find_nsectionsdat(path)
 xxx = read_settings(path,nsectionsdat)

 
; if (itp lt xxx.min) then begin
;  print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
;  stop
; endif
 
 RR = read_major_radius(path)
;#; getting the files I need
;#; I need B, v
 if ( (pdbr eq !null) eq 1) then begin
 restore_ifile = path+'dat/ibprofiles.sav'
 ex_file = file_test(restore_ifile)
 if (ex_file ne 1) then begin
    print,'need to create the fields file ibprofiles.sav',d3
    call_procedure,'isavemodeprofiles',path,'0',d3
 endif
 restore,restore_ifile
 endif

 if ( (pdvr eq !null) eq 1) then begin
 restore_ifile = path+'dat/ivprofiles.sav'
 ex_file = file_test(restore_ifile)
 if (ex_file ne 1) then begin
    print,'need to create the fields file ivprofiles.sav',d3
    call_procedure,'isavemodeprofiles',path,'0',d3
 endif
 restore,restore_ifile
 endif

 if ( (pjr eq !null) eq 1) then begin
;#; I need J
 restore_ifile = path+'dat/ijprofiles.sav'
 ex_file = file_test(restore_ifile)
 if (ex_file ne 1) then begin
    print,'need to create the fields file ijprofiles.sav',d3
    call_procedure,'isavejprofiles',path,'0',d3
 endif
 restore,restore_ifile
 endif

 elr = reform(pdbt)*0.d
 tin = 0
 tfin = n_elements(pdbr[*,0,0])-1

;#; first: interpolation of br,vr on the x2 mesh
; if ( (ddbr eq !null) eq 1) then begin
 ddbr =  elr * dcomplex(0.d,0.d)
 ddvr =  elr * dcomplex(0.d,0.d)
 ddjt =  elr * dcomplex(0.d,0.d)
 ddjz =  elr * dcomplex(0.d,0.d)
 help,square
 print,nanf
 for q=tin, tfin do begin
  for j=0, n_elements(elr[0,0,*])-1 do begin
        dbr = interpol(reform(pdbr[q,*,j]), pror1, pror2 )
        dvr = interpol(reform(pdvr[q,*,j]), pror1, pror2 )
        djt = interpol(reform( pjt[q,*,j]), pror1, pror2 )
        djz = interpol(reform( pjz[q,*,j]), pror1, pror2 )
        mn = mnum(j,nz,nanf,d3)
        if (mn[0] eq 0) then begin
         dbr[0] = -dbr[1]
         dvr[0] = -dvr[1]
         djt[0] = -djt[1]
         djz[0] = djz[1]
        endif else begin
         if (mn[1] eq 1) then begin
          dbr[0] = dbr[1]
          dvr[0] = dvr[1]
          djt[0] = djt[1]
          djz[0] = -djz[1]
         endif else begin
          dbr[0] = -dbr[1]
          dvr[0] = -dvr[1]
          djt[0] = -djt[1]
          djz[0] = djz[1]
         endelse
        endelse
    ddbr[q,*,j] = dbr
    ddvr[q,*,j] = dvr
    ddjt[q,*,j] = djt
    ddjz[q,*,j] = djz
  endfor
  endfor
  print,'br, vr interpolated'
;  endif

; tin=0
; tfin = 210
 jin = 70
 jfin = 80
 eta = 1.e-6

 if (equ eq 1) then begin
;#; compute only axysymmetric electric field

  elt = elr
  elz = elr
  print,'computing only equilibrium'
  for q=tin, tfin do begin
  if (q mod 10 eq 0) then begin
   print,'computing ITP= ',q
  endif
   btvz =  nonlinearconvo(path,reform(pdbt[q,*,*]),reform(pdvz[q,*,*]),0,d3)
   bzvt =  nonlinearconvo(path,reform(pdbz[q,*,*]),reform(pdvt[q,*,*]),0,d3)
   bzvr =  nonlinearconvo(path,reform(pdbz[q,*,*]),reform(ddvr[q,*,*]),0,d3)
   brvz =  nonlinearconvo(path,reform(ddbr[q,*,*]),reform(pdvz[q,*,*]),0,d3)
   brvt =  nonlinearconvo(path,reform(ddbr[q,*,*]),reform(pdvt[q,*,*]),0,d3)
   btvr =  nonlinearconvo(path,reform(pdbt[q,*,*]),reform(ddvr[q,*,*]),0,d3)
   elr[q,*,0] = btvz - bzvt + eta*pjr[q,*,0]
   elt[q,*,0] = bzvr - brvz + eta*ddjt[q,*,0]
   elz[q,*,0] = brvt - btvr + eta*ddjz[q,*,0]
  endfor
  eq_elr = reform(elr[tin:tfin,*,0])
  eq_elt = reform(elt[tin:tfin,*,0])
  eq_elz = reform(elz[tin:tfin,*,0])
  if (mf eq 1) then begin
   out=path+'dat/imf_el_eq.sav'
   strg = 'el'
   call_procedure,'mf_s',pdt,eq_elr,eq_elt,eq_elz,out,strg
  endif else begin
   save_eqfile = path+'dat/el_eq.sav'
   save,filename= save_eqfile,pror1,pror2,eq_elr,eq_elt,eq_elz,tin,tfin
  endelse

 endif else begin

 elt = elr
 elz = elr
 for q=tin, tfin do begin
 print,'ATTENZIONE! calcolo solo uno spettro ristretto!'
;#; equilibrium part
  if (q mod 10 eq 0) then begin
   print,'computing ITP= ',q
  endif
  j=0
  btvz =  nonlinearconvo(path,reform(pdbt[q,*,*]),reform(pdvz[q,*,*]),j,d3)
  bzvt =  nonlinearconvo(path,reform(pdbz[q,*,*]),reform(pdvt[q,*,*]),j,d3)
  bzvr =  nonlinearconvo(path,reform(pdbz[q,*,*]),reform(ddvr[q,*,*]),j,d3)
  brvz =  nonlinearconvo(path,reform(ddbr[q,*,*]),reform(pdvz[q,*,*]),j,d3)
  brvt =  nonlinearconvo(path,reform(ddbr[q,*,*]),reform(pdvt[q,*,*]),j,d3)
  btvr =  nonlinearconvo(path,reform(pdbt[q,*,*]),reform(ddvr[q,*,*]),j,d3)
  elr[q,*,j] = btvz - bzvt + eta*pjr[q,*,0]
  elt[q,*,j] = bzvr - brvz + eta*pjt[q,*,0]
  elz[q,*,j] = brvt - btvr + eta*pjz[q,*,0]
  for j=jin, jfin do begin
  btvz =  nonlinearconvo(path,reform(pdbt[q,*,*]),reform(pdvz[q,*,*]),j,d3)
  bzvt =  nonlinearconvo(path,reform(pdbz[q,*,*]),reform(pdvt[q,*,*]),j,d3)
  bzvr =  nonlinearconvo(path,reform(pdbz[q,*,*]),reform(ddvr[q,*,*]),j,d3)
  brvz =  nonlinearconvo(path,reform(ddbr[q,*,*]),reform(pdvz[q,*,*]),j,d3)
  brvt =  nonlinearconvo(path,reform(ddbr[q,*,*]),reform(pdvt[q,*,*]),j,d3)
  btvr =  nonlinearconvo(path,reform(pdbt[q,*,*]),reform(ddvr[q,*,*]),j,d3)
;  if (q eq 200 and j eq 73) then begin
;   plot,btvz,yr=[-0.01,0.01]
;   oplot,bzvt,col=50,psym=symcat(16)
;   oplot,bzvr,col=100,psym=symcat(16)
;   oplot,brvz,col=150,psym=symcat(16)
;   oplot,brvt,col=200,psym=symcat(16)
;   oplot,btvr,col=250,psym=symcat(16)
;  endif
   elr[q,*,j] = btvz - bzvt + eta*pjr[q,*,j]
   elt[q,*,j] = bzvr - brvz + eta*ddjt[q,*,j]
   elz[q,*,j] = brvt - btvr + eta*ddjz[q,*,j]

 endfor
 endfor
 save_eqfile = path+'dat/el.sav'
 save,filename= save_eqfile,pror1,pror2,elr,elt,elz,tin,tfin,jin,jfi
 endelse

; restore_ifile = path+'dat/ijprofiles.sav'
; spawn,'rm -f '+restore_ifile
; restore_ifile = path+'dat/ibprofiles.sav'
; spawn,'rm -f '+restore_ifile
; restore_ifile = path+'dat/ivprofiles.sav'
; spawn,'rm -f '+restore_ifile

 endprogram: print,'sono alla fine, non serviva ri-calcolare il campo elettrico'

end
