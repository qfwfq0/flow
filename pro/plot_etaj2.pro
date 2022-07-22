pro plot_etaj2, path, d3, tok, nth, nzz, itp0, itp1, hel, deltaitp, mymin, mymax


 jparfile = path+'dat/jpar_j2_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'_itp'+strtrim(string(itp0,format='(i0)'))+strtrim(string(itp1,format='(i0)'))+'.sav'
 result = file_test(jparfile)
 if (result ne 1) then begin
  call_procedure, 'isave_jpar_j2', path, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
 endif else begin
  restore,jparfile
 endelse

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
  print,'dir',dir
  call_procedure,'check_dir',dir 
 endfor

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
  fexis=file_test(file)
  if (fexis ne 1) then begin
  call_procedure, 'time_chihel',path,d3,tok,hel,itp0,itp1,deltaitp
  endif
  restore,file
 endif

 
   rchi = reform(rchi)
 help,nsaved
 
 if (tok eq 1) then begin
;   image_file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/flowxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
;   letter="161B;"
;   greektheta='!9' + String(letter) + '!X'
;   letter2="160B;"
;   greekpi='!9' + String(letter2) + '!X'
;   letter="143B;"
;   greekchi='!9' + String(letter) + '!X'
;   yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
;   ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]
;
;   !p.psym=-3
;;;#; interpolo br e vr su mesh x2
;   help,vr
;   help,br
;    dbr = bt*0.d
;    dvr = bt*0.d
;    for l=0, n_elements(br[0,0,*])-1 do begin
;    for j=0, n_elements(br[0,*,0])-1 do begin
;     dbr[*,j,l] = interpol(reform(br[*,j,l]), pror1, pror2 )
;     dvr[*,j,l] = interpol(reform(vr[*,j,l]), pror1, pror2 )
;    endfor
;    endfor
;
;;#; scelgo il z
;    zzz = 0
;;#; scrivo nuova mesh xx,yy
;     nn = 51
;     xx = make_array(nn,/double,value=0.d)
;     yy = make_array(nn,/double,value=0.d)
;     vx = make_array(nn,nn,/double,value=0.d)
;     vy = make_array(nn,nn,/double,value=0.d)
;     xx = findgen(nn)/(nn-1)*2.-1.
;     yy = findgen(nn)/(nn-1)*2.-1.
;     dth = vecth(5) - vecth(4)
;     dr = pror2(5) - pror2(4)
;     for i=0,nn-1 do begin
;      for k=0,nn-1 do begin
;       if( (xx(i)^2.+yy(k)^2.) le 1. ) then begin
;;#; calcolo r,th corrispondenti ai miei punti in x,y
;        rr = sqrt(xx(i)^2.+yy(k)^2.)
;        th = atan(yy(k),xx(i))
;        if (th lt 0.d) then th = th + 2.*!pi
;;#; calcolo gli indici frazionari per interpolare
;        loc_th = value_locate(vecth,th)
;        idx_th = loc_th + (th-vecth(loc_th))/dth
;        loc_r = value_locate(pror2,rr)
;        idx_r = loc_r + (rr-pror2(loc_r))/dr
;;        print,idx_th,idx_r
;;#; interpolo 
;        tvr = interpolate(dvr[*,*,zzz],idx_r,idx_th,cubic=-0.5)
;        tvt = interpolate(vt[*,*,zzz],idx_r,idx_th,cubic=-0.5)
;        vx(i,k) = tvr*cos(th) - tvt*sin(th)
;        vy(i,k) = tvr*sin(th) + tvt*cos(th)
;       endif else begin
;        vx(i,k) = 0.d
;       endelse
;      endfor
;     endfor
;     
;;#; find separatrix
;   q=0
;   aa = min(rchi[*,*],idx) ;#; O-point
;   bb = array_indices(rchi[*,*],idx)
;   cc = max(rchi[bb[0],*]) ;#; separatrix!
;   print,'sep_lev=',cc
;   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
;
;   deltachi = max(rchi)-min(rchi)
;
;   nisla=4
;   dd=findgen(nisla)/(nisla-1)*abs((cc-aa))*0.75+aa
;   print,'dd',dd
;;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;;   print,'ee',ee
;   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    xrange=[-1.,1.]
;    yrange = xrange
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;
;    xyouts,0.6,0.9,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7
;    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
;   phd_graphcloser,image_file
;
;
;   loadct,13
;   jr_file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/jrxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
;   phd_graphopener,filename=jr_file,xsize=12.,ysize=12.
;    xrange=[-1.,1.]
;    yrange = xrange
;    contour,reform(jr[*,*,0]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!Dr!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
;    contour,reform(jr[*,*,0]),pror2#cos(vecth),pror2#sin(vecth),nlev=200,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;   phd_graphcloser,jr_file
;   jt_file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/jtxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
;   phd_graphopener,filename=jt_file,xsize=12.,ysize=12.
;    xrange=[-1.,1.]
;    yrange = xrange
;    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!D'+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
;    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=200,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;   phd_graphcloser,jt_file
;   jz_file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/jzxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
;   phd_graphopener,filename=jz_file,xsize=12.,ysize=12.
;    xrange=[-1.,1.]
;    yrange = xrange
;    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!Dz!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
;    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=200,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;   phd_graphcloser,jz_file
;
;
;infiles = [jr_file,jt_file,jz_file]
;infiles2 = strjoin(infiles, ' ')
;outfile = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/currentmyminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
;spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
;for i=0,n_elements(infiles) -1 do begin
; spawn, "rm -f "+infiles(i)
; print,infiles(i)
;endfor

endif else begin ;#; begin RFP part

  print,'nsaved= ',nsaved
  print,'ITP',iitp
 
 for q=0,n_elements(iitp)-1 do begin
   print,'avanti sempre, q=',q,nsaved,iitp(q)
   time_file = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/etaj2_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
   restore,time_file
   image_file = path+'dat/itp/'+strtrim(string(iitp(q),'(i0)'))+'/etaj2_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   letter="161B;"
   greektheta='!9' + String(letter) + '!X'
   letter2="160B;"
   greekpi='!9' + String(letter2) + '!X'
   letter="143B;"
   greekchi='!9' + String(letter) + '!X'
   letter="150B;"
   greeketa='!9' + String(letter) + '!X'
   yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
   ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

   !p.psym=-3
     
;#; find separatrix
   aa = max(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   ee0 = min(rchi[bb[0],*],eeidx)
   ee = max(rchi[*,eeidx])
   print,'sep_lev=',ee
   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
   deltachi = max(rchi)-min(rchi)

   nisla=5
   dd=findgen(nisla)/(nisla-1)*abs((ee-aa))+ee
   print,'dd',dd
;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;   print,'ee',ee
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title=greekchi+'!Dhel!N   h='+strtrim(string(hel,'(i0)')),xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=18,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot

    xyouts,0.6,0.9,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7
   phd_graphcloser,image_file


;#; create resistivity vector
  eta_vec = make_array(n_elements(pror2),n_elements(ajsq[0,*,0]),/double,value=0.d)
  for iq=0,n_elements(eta_vec[0,*])-1 do begin
   eta_vec[*,iq] = dissipation.eta*(1+dissipation.alet*pror2^(dissipation.beet))^dissipation.gaet
  endfor
;   loadct,,,
   j2_file = path+'dat/itp/'+strtrim(string(iitp(q),'(i0)'))+'/etaj2_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=j2_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(ajsq[*,*,0])*eta_vec,pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title=greeketa+'J!E2!N   h='+strtrim(string(hel,'(i0)')),xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(ajsq[*,*,0])*eta_vec,pror2#cos(vecth),pror2#sin(vecth),nlev=25,/fill,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=184,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
    xyouts,0.6,0.96,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7
   phd_graphcloser,j2_file

   jpar_file = path+'dat/itp/'+strtrim(string(iitp(q),'(i0)'))+'/jpar_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jpar_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(ajpar[*,*,0]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title='J!Epar!N   h='+strtrim(string(hel,'(i0)')),xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(ajpar[*,*,0]),pror2#cos(vecth),pror2#sin(vecth),nlev=25,/fill,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=184,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
    xyouts,0.6,0.96,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7
   phd_graphcloser,jpar_file

  loadct,0,file='~/idl/colors.tbl'
   jpar2_file = path+'dat/itp/'+strtrim(string(iitp(q),'(i0)'))+'/jpar2_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jpar2_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    xrange=[-0.75,0.75]
    yrange = xrange
;#; specifico per un caso paper riconnessione 2020 / simulazione vario_dissipazione/nmp7_4.5e-3/
;    levels=(findgen(25)/24-0.5)*16.5+4.
;    if (q eq 0) then begin
;     print,'levels',levels
;    endif
;    zmax=464
;#; specifico per un caso paper riconnessione 2020 / simulazione 2D h9 s100m1000 qinfty momsour 0
    levels=(findgen(25)/23-0.5)*12.5+4.
;    if (q eq 0) then begin
;     print,'levels',levels
;    endif
;    zmax=66

    qqzmax = max(ajpar,idx)
    qzmax=array_indices(ajpar,idx)
    zmax = qzmax[2]
;    zmax=vecz(qzmax[2])/4.
;    contour,reform(ajpar[*,*,zmax]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title='z/R!D0!N='+strtrim(string(zmax,'(f7.3)')),xtitle='x/a',ytitle='y/a',/nodata,/iso,chars=1.7
;    zmax=
    mxjpa=max(ajpar,idx)
    idjpa=array_indices(ajpar,idx)
    print,'ITP=',iitp(q),mxjpa,idjpa(2)
    contour,reform(ajpar[*,*,zmax]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title='',xtitle='x/a',ytitle='y/a',/nodata,/iso,chars=1.7,position=[0.21,0.18,0.95,0.95]
;    contour,reform(ajpar[*,*,zmax]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title='J!Dpar!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,chars=2.
   loadct,33,/silent
    contour,reform(ajpar[*,*,zmax]),pror2#cos(vecth),pror2#sin(vecth),lev=levels,/fill,/overplot
;    contour,reform(ajpar[*,*,zmax]),pror2#cos(vecth),pror2#sin(vecth),nlev=20,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;loadct,7
;#; specifico per un caso paper riconnessione 2020 / simulazione 2D h9 s100m1000 qinfty momsour 0
;   rchi=shift(rchi,0,-43) 
;    contour,reform(rchi[*,*]),pror1#cos(vectth-1.0716),pror1#sin(vectth-1.0716),lev=ee[0],col=100,thick=2.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=100,thick=2.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd[3],col=25,thick=1.,/overplot
;;;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=25,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=25,thick=1.,/overplot
;   loadct,33,/silent
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    xyouts,0.75,0.12,/norm, 'ITP='+strtrim(string(itp(q)),'(i0)'),chars=0.7
;    xyouts,0.5,0.81,/norm, 'ITP='+strtrim(string(itp(q)),'(i0)'),chars=1.
loadct,0,/silent
    xyouts,0.58,0.95,/norm, 'J!Dpar!N',chars=1.

loadct,33,/silent
   phd_graphcloser,jpar2_file
  loadct,0,file='~/idl/colors.tbl'

    jparprof_file = path+'dat/itp/'+strtrim(string(iitp(q),'(i0)'))+'/jpar_profile_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
    phd_graphopener,filename=jparprof_file,xsize=12.,ysize=12.
    levels=(findgen(25)/24-0.5)*16.5+4.
    if (q eq 0) then begin
     print,'levels',levels
    endif
;#; specifico per un caso paper riconnessione 2020 / simulazione vario_dissipazione/nmp7_4.5e-3/
;    zmax=464
;    thmax=30 ;#; specifico per un caso
   thmax=5 ;#; specifico per un caso
;    thmax=qzmax[1]
    thmax2=thmax+n_elements(ajpar[0,*,0])/2 ;#; specifico per un caso
    xrange2=[-1.,1.]
    yrange2=[-2.,10.]
    radius=[-reverse(pror2[2:-1]),pror2]
    plot,[0],xrange=xrange2,yrange=yrange2,xstyle=1,ystyle=1,xtitle='r/a',ytitle='J!Dpar!N',title='',/nodata,col=0,chars=1.7,position=[0.21,0.18,0.95,0.95]
    oplot,pror2,ajpar[*,thmax,zmax],thick=6.,col=40
    oplot,-reverse(pror2),reverse(ajpar[*,thmax2,zmax]),thick=6.,col=40
;    xyouts,0.5,0.32,/norm, 'ITP='+strtrim(string(iitp(q)),'(i0)'),chars=1.
    xyouts,0.5,0.37,/norm, 'J!Dpar!N tot',chars=1.3,col=40
;#; per plottare la parte assialsimmetrica
;    xyouts,0.5,0.30,/norm, 'J!Dpar!N m=0,n=0',chars=1.3,col=0
;    j00_file = path+'dat/itp/'+strtrim(string(iitp(q),format='(i0)'))+'/etaj2_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax00.sav'
;    restore,j00_file
;    oplot,pror2,ajpar[*,0,0],thick=3.,col=0,linestyle=0
;    oplot,-reverse(pror2),reverse(ajpar[*,0,0]),thick=3.,col=0,linestyle=0
;    oplot,pror2,deriv(pror2,ajpar[*,0,0])/5.,thick=3.,col=250,linestyle=0

    
   phd_graphcloser,jparprof_file

  loadct,0,file='~/idl/colors.tbl'
   
   jparz_file = path+'dat/itp/'+strtrim(string(iitp(q),'(i0)'))+'/jparz_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jparz_file,xsize=12.,ysize=12.
    xrange=[0.,1.]
    yrange = [0.,8.*!pi]
    qqthmax = max(ajpar,idx)
    qthmax=array_indices(ajpar,idx)
    thmax = qthmax[1]
    contour,reform(ajpar[*,thmax,*]),pror2,vecz,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title='J!Epar!N   h='+strtrim(string(hel,'(i0)'))+'  ITP='+strtrim(string(iitp(q)),'(i0)'),xtitle='x/a',ytitle='y/a',/nodata
    contour,reform(ajpar[*,thmax,*]),pror2,vecz,lev=levels,/fill,/overplot
    xyouts,0.6,0.96,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7

   phd_graphcloser,jparz_file

;   jt_file = path+'dat/itp/'+strtrim(string(iitp0,'(i0)'))+'/jtxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
;   phd_graphopener,filename=jt_file,xsize=12.,ysize=12.
;    xrange=[-1.,1.]
;    yrange = xrange
;    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title='J!D'+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
;    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=25,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    xyouts,0.6,0.9,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7
;   phd_graphcloser,jt_file
;   jz_file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/jzxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
;   phd_graphopener,filename=jz_file,xsize=12.,ysize=12.
;    xrange=[-1.,1.]
;    yrange = xrange
;    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=25,title='J!Dz!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
;    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=25,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    xyouts,0.6,0.9,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7
;   phd_graphcloser,jz_file
;

 infiles = [j2_file,jpar_file]
 infiles2 = strjoin(infiles, ' ')
 outfile = path+'dat/itp/'+strtrim(string(iitp(q),'(i0)'))+'/j2jpar_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
 spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
 for i=0,n_elements(infiles) -1 do begin
  spawn, "rm -f "+infiles(i)
  print,infiles(i)
 endfor
;print,'fatto q= ',q,nsaved
endfor ;#; end of cycle on nsaved
endelse
end
