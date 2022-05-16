pro flow_xy,path,d3,tok,hel,itp0,itp1,deltaitp

;#; check chosen times
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
 itp= indgen(nsaved)*deltaitp + itp0 
 print,'itp',itp
 for q=0,n_elements(itp)-1 do begin
  dir =  path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/'
  print,'dir',dir
  call_procedure,'check_dir',dir 
 endfor
   
   ;path = './'
 if (rchi eq !null) then begin
  file=path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/chihel.sav'
  fexis=file_test(file)
  if (fexis ne 1) then begin
  call_procedure, 'time_chihel',path,d3,tok,hel,itp0,itp1,deltaitp
  endif
  restore,file
 endif

 if (vr eq !null) then begin
  file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/antifft_v_102x4.sav'
  fexis=file_test(file)
  if (fexis ne 1) then begin
   if (tok) then begin
    str='v'
    nth = 102
    nzz = 4
    call_procedure, 'antifft_reconst_tok', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp
   endif
  endif
  restore,file
  file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/antifft_b_102x4.sav'
  fexis=file_test(file)
  if (fexis ne 1) then begin
   if (tok) then begin
    str='b'
    nth = 102
    nzz = 4
    call_procedure, 'antifft_reconst_tok', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp
   endif
  endif
  restore,file
 endif

   rchi = reform(rchi)
 
   image_file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/flowxy.eps'
   letter="161B;"
   greektheta='!9' + String(letter) + '!X'
   letter2="160B;"
   greekpi='!9' + String(letter2) + '!X'
   letter="143B;"
   greekchi='!9' + String(letter) + '!X'
   yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
   ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

   !p.psym=-3
;;#; interpolo br e vr su mesh x2
   help,vr
   help,br
    dbr = bt*0.d
    dvr = bt*0.d
    for l=0, n_elements(br[0,0,*])-1 do begin
    for j=0, n_elements(br[0,*,0])-1 do begin
     dbr[*,j,l] = interpol(reform(br[*,j,l]), pror1, pror2 )
     dvr[*,j,l] = interpol(reform(vr[*,j,l]), pror1, pror2 )
    endfor
    endfor

;#; scelgo il z
    zzz = 0
;#; scrivo nuova mesh xx,yy
     nn = 51
     xx = make_array(nn,/double,value=0.d)
     yy = make_array(nn,/double,value=0.d)
     vx = make_array(nn,nn,/double,value=0.d)
     vy = make_array(nn,nn,/double,value=0.d)
     xx = findgen(nn)/(nn-1)*2.-1.
     yy = findgen(nn)/(nn-1)*2.-1.
     dth = vecth(5) - vecth(4)
     dr = pror2(5) - pror2(4)
     for i=0,nn-1 do begin
      for k=0,nn-1 do begin
       if( (xx(i)^2.+yy(k)^2.) le 1. ) then begin
;#; calcolo r,th corrispondenti ai miei punti in x,y
        rr = sqrt(xx(i)^2.+yy(k)^2.)
        th = atan(yy(k),xx(i))
        if (th lt 0.d) then th = th + 2.*!pi
;#; calcolo gli indici frazionari per interpolare
        loc_th = value_locate(vecth,th)
        idx_th = loc_th + (th-vecth(loc_th))/dth
        loc_r = value_locate(pror2,rr)
        idx_r = loc_r + (rr-pror2(loc_r))/dr
;        print,idx_th,idx_r
;#; interpolo 
        tvr = interpolate(dvr[*,*,zzz],idx_r,idx_th,cubic=-0.5)
        tvt = interpolate(vt[*,*,zzz],idx_r,idx_th,cubic=-0.5)
        vx(i,k) = tvr*cos(th) - tvt*sin(th)
        vy(i,k) = tvr*sin(th) + tvt*cos(th)
       endif else begin
        vx(i,k) = 0.d
       endelse
      endfor
     endfor
     
;#; find separatrix
   q=0
   aa = min(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   cc = max(rchi[bb[0],*]) ;#; separatrix!
   print,'sep_lev=',cc
   print,'Opoint_lev=',aa;,min(cc[0],aa[0])

   deltachi = max(rchi)-min(rchi)

   nisla=4
   dd=findgen(nisla)/(nisla-1)*abs((cc-aa))*0.75+aa
   print,'dd',dd
;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;   print,'ee',ee
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot

    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5


;    maxvz=max(vz[*,*,zzz])
;;    contour,vz[*,*,zzz]/maxvz*deltachi+min(rchi),pror2,vecth,nlev=40,c_col=20,/overplot,/data
;;    contour,vz[*,*,zzz]/maxvz*deltachi+min(rchi),pror2,vecth,nlev=40,/fill,/overplot,/data
;
;;    maxvr=max(vr[*,*,zzz])
;;    minvr=min(vr[*,*,zzz])
;;    cacca=[abs(minvr[0]),abs(maxvr[0])]
;;    maxvr=max(cacca)
;;    newvr=(vr[*,*,zzz]-minvr)
;;    maxnewvr=max(newvr)
;;    newvr= newvr/maxnewvr*deltachi+min(rchi)
;;    contour,newvr,pror1,vecth,nlev=40,/fill,/overplot,/data
;
;    maxvt=max(vt[*,*,zzz])
;    contour,vt[*,*,zzz]/maxvt*deltachi+min(rchi),pror2,vecth,nlev=40,/fill,/overplot,/data
;;    contour,vt[*,*,zzz]/maxvt*deltachi+min(rchi),pror2,vecth,nlev=40,/overplot,/data
;
;;    maxbr=max(br[*,*,zzz])
;;    contour,br[*,*,zzz]/maxbr*deltachi+min(rchi),pror1,vecth,nlev=40,/overplot,/data
;
;    contour,reform(rchi[*,*]),pror1,vectth,lev=cc[0],col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1,vectth,levels=dd,col=253,thick=1.,/overplot

   phd_graphcloser,image_file
end
