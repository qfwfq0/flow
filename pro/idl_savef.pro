;function janz,m,n,my,mm,nanf,nz
;    j = 0
;    if (m eq 0 and n eq 0) then return,j
;    for mi=0,m-1 do begin
;        j = j + nz[mi]
;    endfor
;    j =  j + 1 + n - nanf[m]
;    return,j
;end

pro idl_savef,path,d3
;#; this routine computes the reversal parameter and the pinch parameter in time

  RR = read_major_radius(path)

;  restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
;  idl_file = path+'dat/ib00profiles.sav'
  mf = 1
  idl_file = path+'dat/b_eq.sav'
  call_procedure,'check_sections',path,mf,d3
  exist_file = file_test(idl_file)
  if (exist_file ne 1) then begin
   print,'need to create the b_eq.sav file!'
;   call_procedure,'isaveeq',path
   exist2 = file_test(path+'dat/imf_bprofiles.sav')
   if (exist2 ne 1) then begin
    print,'fare qualcosa!!!'
;!    call_procedure,'ireadeq',path
   endif else begin
    restore,path+'dat/imf_bprofiles.sav'
    eqt = tt
    eqbr = reform(br[*,*,0,0]) * sin(reform(br[*,*,0,1])) 
    eqbt = reform(bt[*,*,0,0]) * cos(reform(bt[*,*,0,1])) 
    eqbz = reform(bz[*,*,0,0]) * cos(reform(bz[*,*,0,1])) 
    lx = n_elements(pror1) - 1
    save,filename=path+'dat/b_eq.sav',eqt,lx,eqbr,eqbt,eqbz
   endelse
  endif
  restore,idl_file

  nr = n_elements(eqt)
  lx = 100
  pror1 = findgen(lx+1) / lx
  pror2 = findgen(lx+2) / lx - 1./(2.*lx)
  phitor = make_array(n_elements(eqt),n_elements(pror2)-1,/double,value=0.) 
  dr = 1.d / (n_elements(pror2)-2)

  for p=0,n_elements(phitor(*,0))-1 do begin
  for q=1,n_elements(phitor(0,*))-1 do begin
   phitor(p,q) = phitor(p,q-1) + dr * pror2(q) * eqbz(p,q)
  endfor
  endfor

  phitor = 2. *!Dpi * phitor
  
;#; reversal parameter
  ff = make_array(n_elements(eqt),/double,value=0.) 
  theta = make_array(n_elements(eqt),/double,value=0.) 
  
  ff = eqbz[*,n_elements(eqbz[0,*])-1] / ( (phitor(*,n_elements(phitor[0,*])-1)) / (!Dpi))
  theta = eqbt[*,n_elements(eqbz[0,*])-1] / ( (phitor(*,n_elements(phitor[0,*])-1)) / (!Dpi))

  sav_file = path+'dat/f_th.sav'
  save,filename=sav_file,eqt,pror1,ff,theta,phitor

;#; plots
  filename = path+'dat/f_th.eps'
  xsize=20
  ysize=12
  phd_graphopener,filename=filename,xsize=xsize,ysize=ysize
  cs=2.
  loadct,0,file='~/idl/colors.tbl'
   
  letter="164B;"
  greektau='!9' + String(letter) + '!X'
  letter="121B;"
  bigtheta='!9' + String(letter) + '!X'
  m0 = 0.12
  m1 = 0.95
  m10 = 0.12
  m11= 0.49
  m12 = 0.57
  m13 = 0.95

  !p.charsize=1.2
   !p.multi=[0,1,2,0,0]
   xr=[min(eqt),max(eqt)]
   yfr=[-0.4,0.1]
   plot,eqt,ff,/nodata,xrange=xr,yrange=yfr,xstyle=1,position=[m0,m12,m1,m13],ytitle='F'
   oplot,eqt,ff,thick=6.,col=200

   !p.multi=[1,1,2,0,0]
   ythr=[1.3,1.7]
   plot,eqt,theta,/nodata,xrange=xr,yrange=ythr,xstyle=1,position=[m0,m10,m1,m11],xtitle='time ('+greektau+'!DA!N)',ytitle=bigtheta
   oplot,eqt,theta,thick=6.,col=10
  phd_graphcloser,filename

end
