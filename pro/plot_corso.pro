;#; Marco, plot per corso nonlinear MHD marzo 2021
pro plot_corso, path, d3, tok, nth, nzz, itp0, itp1, hel, deltaitp, mymin, mymax


;#; check chosen times
 dissipation = read_dissipation(path)
 print,'dissipation',dissipation.eta,dissipation.nu
 if (dissipation.eta*dissipation.nu le 0.d) then begin 
  print,'si è letta male la dissipazione'
  stop
 endif
  
; print,dissipation
 print,'AAA, tmin,tmax',itp0,itp1
 restore,path+'dat/imf_bprofiles.sav'
 if (itp1 eq -1) then begin 
  itp1=n_elements(br[*,0,0,0])-1
  itp0=n_elements(br[*,0,0,0])-1
  deltaitp=1
 endif
 if (itp0 eq itp1) then begin
;  outfile = path+'dat/itp/'+strtrim(string(itp0,format='(i0)'))+'/chihel.sav'
  print,'need to create the directory for helical flux function file for itp=',itp0
  dir =  path+'dat/itp/'+strtrim(string(itp0,format='(i0)'))+'/'
  print,'dir',dir
  call_procedure,'check_dir',dir 
 endif
 itp0 = fix(itp0,type=3)
 itp1 = fix(itp1,type=3)
 nitp = itp1 - itp0 + 1
 nsaved = ceil(float(nitp)/deltaitp)
 iitp= indgen(nsaved)*deltaitp + itp0 
 print,'itpAAA',iitp
 for q=0,n_elements(iitp)-1 do begin
  dir =  path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/'
;  print,'dir',dir
  call_procedure,'check_dir',dir 
 endfor

;#; mi servono fattore di sicurezza, jpar, helical flux function, velocityfield
 jparfile = path+'dat/jpar_j2_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'_itp'+strtrim(string(itp0,format='(i0)'))+strtrim(string(itp1,format='(i0)'))+'.sav'
 result = file_test(jparfile)
 if (result ne 1) then begin
  call_procedure, 'isave_jpar_j2', path, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
 endif else begin
  restore,jparfile
 endelse



 result = file_test(path+'dat/spectrum.sav')
 if (result ne 1) then begin
  call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
  print,'nanf[1]',nanf[1]
 endif else begin
  restore,path+'dat/spectrum.sav'
  nzcon=nz
  nzcon[0]=nz[0]+1
 endelse

 if (d3 eq 0) then begin
  hel = nanf(1)
 endif else begin
  print,'cosa vale hel?',hel
  if (hel eq !null) then stop
 endelse

   ;path = './'
 if (rchi eq !null) then begin
  file=path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/chihel.sav'
;  fexis=file_test(file)
;  if (fexis ne 1) then begin
  call_procedure, 'time_chihel',path,d3,tok,hel,itp0,itp1,deltaitp
;  endif
  restore,file
 endif

 
 rchi = reform(rchi)
 help,nsaved
 if (vr eq !null) then begin
  file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/antifft_v_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,'(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
;  fexis=file_test(file)
;  if (fexis ne 1) then begin
   if (tok eq 1) then begin
    str='v'
    nth = nth
    nzz = nzz
    call_procedure, 'antifft_reconst_tok', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp, mymin, mymax
   endif else begin
    str='v'
    nth = nth
    nzz = nzz
    call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
   endelse
;  endif
  restore,file
 endif
  
;#; Marco, chech q profiles
 qprof_file = path+'/dat/qprof.sav'
 qfile_ex = file_test(qprof_file)
 if (qfile_ex ne 1) then begin
  call_procedure, 'qprof',path
 endif else begin
  restore,qprof_file
 endelse

!p.charsize=2.
 itp0 = fix(itp0,type=3)
 itp1 = fix(itp1,type=3)
 nitp = itp1 - itp0 + 1
 nsaved = ceil(float(nitp)/deltaitp)
 fiitp= indgen(nsaved)*deltaitp + itp0 
 print,'itpAAA2',fiitp
help,fiitp
 for q=0,n_elements(fiitp)-1 do begin
 print,'making ',fiitp(q),' time step'
;#; restore Jpar
 jparfile = path+'dat/itp/'+strtrim(string(fiitp(q),'(i0)'))+'/etaj2_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
 restore,jparfile

;#; Immagine
 image_file = path+'dat/itp/'+strtrim(string(fiitp(q),'(i0)'))+'/corso.eps'
 phd_graphopener,filename=image_file,xsize=17.,ysize=6.
 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter2="160B;"
 greekpi='!9' + String(letter2) + '!X'
 letter="143B;"
 greekchi='!9' + String(letter) + '!X'
 yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
 ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

 m1=0.07 & m2 =0.25 & m3=0.32 & m4=0.5 & m5 = 0.57 & m6 = 0.75 & m7=0.82 & m8=0.99
 m9 = 0.17 & m10 = 0.88
; xyouts,0.34,0.955,'time = '+strtrim(string(ttt,format='(i0)'))+' '+greektau+'!DA!N, ITP='+strtrim(string(itp,'(i0)')),/norm,chars=1.2

 lcol=210
 thkk = 3.
;#; plot qprofile
print,'CHECK',q,fiitp(q)
 !p.multi=[4,4,1,0,0]
 yrq = [-0.1,5.]
 plot,[0.],xrange=xr,yrange=yrq,xtitle='r/a',ytitle='q(r)',title='q(r)',position=[m1,m9,m2,m10],ystyle=1,/nodata,xstyle=1
 oplot,pror2,qprof[fiitp(q),*],col=lcol,psym=-3,thick=thkk
 oplot,pror2,shear_q[fiitp(q),*],col=lcol-200,psym=-3,thick=thkk,linestyle=2
 oplot,pror2,qprof[fiitp(0),*],col=0,psym=-3,thick=0.5,linestyle=2
 xyouts,0.08,0.81,/norm, 'q (orange)',chars=0.7,color=lcol
 xyouts,0.08,0.77,/norm, 'shear_q (blue)',chars=0.7,color=lcol-200
 bb = min(where(qprof[fiitp(q),*] gt 1.))
 r2 = pror2(bb)
 plots,r2,yrq[0],col=1,syms=0.9
; plots,r2,yrq[1],col=1,linestyle=2,thick=1.,syms=1.

 xyouts,0.48,0.95,/norm, 't='+strtrim(string(fiitp(q),'(i0)')),chars=0.9,color=0

;#; plot Jpar
 letter="275B;"
 parallel='!9' + String(letter) + '!X'
 !p.multi=[3,4,1,0,0]
 yrj = [0.,0.7]
 plot,[0.],xrange=xr,yrange=yrj,xtitle='r/a',ytitle='',title='J!D'+parallel+parallel+'!N',position=[m3,m9,m4,m10],ystyle=1,/nodata,xstyle=1
 oplot,pror2,ajpar[*,0,0],col=lcol,psym=-3,thick=thkk
 plots,r2,yrq[0],col=1,syms=0.9
 if (q eq 0) then begin
  jpar0=ajpar[*,0,0]
  oplot,pror2,jpar0,col=0,psym=-3,thick=0.5
 endif else begin
  oplot,pror2,jpar0,col=0,psym=-3,thick=0.5
 endelse

;#; plot helical flux
 chihelfile = path+'dat/itp/'+strtrim(string(fiitp(q),'(i0)'))+'/chihel.sav'
 restore,chihelfile
 !p.multi=[2,4,1,0,0]
 xchi = [0.,1.]
 ychi = [0.,2.*!PI]
 contour,rchi,pror1,vectth,xrange=xchi,yrange=ychi,xtitle='r/a',ytitle=greektheta,title=greekchi+'!Dhel!N',position=[m5,m9,m6,m10],ystyle=1,/nodata,xstyle=1
 contour,rchi,pror1,vectth,nlev=30,c_col=lcol,/overplot
;#; valid for m=1 n=-1
 aaa=max(rchi)
 nlev=5
 bbb=[0.98*aaa,0.985*aaa,0.99*aaa,0.995*aaa,0.997*aaa,0.999*aaa,0.9999*aaa]
 contour,rchi,pror1,vectth,lev=bbb,c_col=lcol+20,/overplot

;#; plot velocity
  vfile = path+'dat/itp/'+strtrim(string(fiitp(q),'(i0)'))+'/antifft_v_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,'(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  restore,vfile
 !p.multi=[1,4,1,0,0]
 contour,vz[*,*,0],pror2,vecth,xrange=xchi,yrange=ychi,xtitle='r/a',ytitle=greektheta,title='v!Dz!N',position=[m7,m9,m8,m10],ystyle=1,/nodata,xstyle=1
 contour,vz[*,*,0],pror2,vecth,nlev=30,/overplot,/fill
; contour,vz[*,*,0],pror2,vecth,nlev=30,c_col=lcol,/overplot
; vr2=vt[*,*,0]*0.
; vr2[0,*] = 0.
; vr2[1:-1,*] = vr[*,*,0]
; velovect,vr2,vt[*:*:4,*:4,0],pror2[*:4],vecth[*:4],length=0.01,thick=0.05,/overplot

 phd_graphcloser,image_file
stop


;#; Immagine
 image_file_th = path+'dat/itp/'+strtrim(string(fiitp(q),'(i0)'))+'/corso_vth.eps'
 phd_graphopener,filename=image_file_th,xsize=17.,ysize=6.
;#; plot qprofile
print,'CHECK',q,fiitp(q)
 !p.multi=[4,4,1,0,0]
 plot,[0.],xrange=xr,yrange=yrq,xtitle='r/a',ytitle='q(r)',title='q(r)',position=[m1,m9,m2,m10],ystyle=1,/nodata,xstyle=1
 oplot,pror2,qprof[fiitp(q),*],col=lcol,psym=-3,thick=thkk
 oplot,pror2,shear_q[fiitp(q),*],col=lcol-200,psym=-3,thick=thkk,linestyle=2
 oplot,pror2,qprof[fiitp(0),*],col=0,psym=-3,thick=0.5,linestyle=2
 xyouts,0.08,0.81,/norm, 'q (orange)',chars=0.7,color=lcol
 xyouts,0.08,0.77,/norm, 'shear_q (blue)',chars=0.7,color=lcol-200
 plots,r2,yrq[0],col=1,syms=0.9
; plots,r2,yrq[1],col=1,linestyle=2,thick=1.,syms=1.

 xyouts,0.48,0.95,/norm, 't='+strtrim(string(fiitp(q),'(i0)')),chars=0.9,color=0

;#; plot Jpar
 !p.multi=[3,4,1,0,0]
 plot,[0.],xrange=xr,yrange=yrj,xtitle='r/a',ytitle='',title='J!D'+parallel+parallel+'!N',position=[m3,m9,m4,m10],ystyle=1,/nodata,xstyle=1
 oplot,pror2,ajpar[*,0,0],col=lcol,psym=-3,thick=thkk
 plots,r2,yrq[0],col=1,syms=0.9
 if (q eq 0) then begin
  oplot,pror2,jpar0,col=0,psym=-3,thick=0.5
 endif else begin
  oplot,pror2,jpar0,col=0,psym=-3,thick=0.5
 endelse

;#; plot helical flux
 !p.multi=[2,4,1,0,0]
 contour,rchi,pror1,vectth,xrange=xchi,yrange=ychi,xtitle='r/a',ytitle=greektheta,title=greekchi+'!Dhel!N',position=[m5,m9,m6,m10],ystyle=1,/nodata,xstyle=1
 contour,rchi,pror1,vectth,nlev=30,c_col=lcol,/overplot
 contour,rchi,pror1,vectth,lev=bbb,c_col=lcol+20,/overplot

;#; plot velocity
 !p.multi=[1,4,1,0,0]
 contour,vt[*,*,0],pror2,vecth,xrange=xchi,yrange=ychi,xtitle='r/a',ytitle=greektheta,title='v!D'+greektheta+'!N',position=[m7,m9,m8,m10],ystyle=1,/nodata,xstyle=1
 contour,vt[*,*,0],pror2,vecth,nlev=30,/overplot,/fill
; contour,vz[*,*,0],pror2,vecth,nlev=30,c_col=lcol,/overplot


 phd_graphcloser,image_file_th


 endfor
end
