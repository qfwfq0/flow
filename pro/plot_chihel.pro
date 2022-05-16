pro plot_chihel, path, d3, tok, dieli,it

  RR = read_major_radius(path)
  if (RR(0) lt 1e-1) then begin
   RR(0)=4.d
  endif

  letter="161B;"
  greektheta='!9' + String(letter) + '!X'
  letter2="160B;"
  greekpi='!9' + String(letter2) + '!X'
  letter="143B;"
  greekchi='!9' + String(letter) + '!X'
  yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
  ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]
  
;#; restore file

  inteli = fix(dieli)
  save_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chi'+strtrim(string(inteli,format='(i0)'))+'.sav'
  restore,save_file
  vectth=th


  phi = 0
  if (tok eq 1) then begin
;#; find separatrix
   q=0
   aa = min(rchi[*,*,phi],idx) ;#; O-point
   bb = array_indices(rchi[*,*,phi],idx)
   cc = max(rchi[bb[0],*,phi]) ;#; separatrix!
;   print,'sep_lev=',cc
;   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
   nisla=8
   dd=findgen(nisla)/(nisla-1)*abs((cc-aa))+aa
;   print,'dd',dd


   image_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel_'+strtrim(string(inteli,format='(i0)'))+'.eps'
;   print,'aaa',image_file
   phd_graphopener,filename=image_file,xsize=15.,ysize=10.
    xrange=[0.,1.]
    yrange=[0.,2.*!pi]
    contour,reform(rchi[*,*]),pror1,vectth,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,/fill,col=0,yticks=5,ytickv=yt0,ytickname=ytpi,title=greekchi+'!Dhel!N '+strtrim(string(it,'(f9.4)')),xtitle='r/a',ytitle=greektheta
    contour,reform(rchi[*,*]),pror1,vectth,lev=cc[0],col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,levels=dd,col=253,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,nlev=50,col=255,thick=0.5,/overplot

   phd_graphcloser,image_file
   image_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel_xy_'+strtrim(string(inteli,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N '+strtrim(string(it,'(f9.4)')),xtitle='x/a',ytitle='y/a',/iso;,position=[0.27,0.21,0.94,0.94]
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=214,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=0,thick=1.,/overplot

   phd_graphcloser,image_file

   image_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel_xy_fill_'+strtrim(string(inteli,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N '+strtrim(string(it,'(f9.4)')),xtitle='x/a',ytitle='y/a',/iso,/fill
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,c_color=255,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=254,thick=3.,/overplot
;   contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=255,thick=1.,/overplot

   phd_graphcloser,image_file

;  endif else begin;#; 2D case
;   print,'code the TOK3D case'
;  endelse;#; 2D case

endif else begin ;#; RFP case
  print,'RFP_case'
   image_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel.eps'
;   print,'image_file',image_file
   letter="161B;"
   greektheta='!9' + String(letter) + '!X'
   letter2="160B;"
   greekpi='!9' + String(letter2) + '!X'
   letter="143B;"
   greekchi='!9' + String(letter) + '!X'
   yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
   ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

;#; altra tecnica per trovare la separatrice guardando il minimo della derivata e poi cercando il massimo della chi a sinistra del minimo della derivata
;   if (bb[0] eq 100) then begin
;    drchi = rchi*0.d
;    for k=0,n_elements(rchi[0,*])-1 do begin
;     drchi[*,k] = deriv(pror1,rchi[*,k])
;    endfor
;     aa=min(drchi,idx)
;     bb=array_indices(drchi,idx)
;   ee0 = min(rchi[bb[0],*],eeidx)
;   ee = max(rchi[*,eeidx])
;   kk=max(rchi[0:bb[0],bb[1]],kidx)
;   endif else begin
;   kk=1.
;   endelse
;
;   print,'sep_lev=',ee,kk
;   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
;   nisla=5
;   dd=findgen(nisla)/(nisla-1)*abs((ee-aa))+ee
;   print,'dd',dd


;#; find separatrix
;   q=0
   aa = max(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   ee0 = min(rchi[bb[0],*],eeidx)
   ee = max(rchi[*,eeidx])
   print,bb
   if (bb[0] eq 100) then begin
    drchi = rchi*0.d
    for k=0,n_elements(rchi[0,*])-1 do begin
     drchi[*,k] = deriv(pror1,rchi[*,k])
    endfor
     aa=min(drchi,idx)
     bb=array_indices(drchi,idx)
   ee0 = min(rchi[bb[0],*],eeidx)
   ropoint = min(rchi[bb[0]:*,bb[1]],idxopoint)
   ee = max(rchi[*,eeidx])
   kk=max(rchi[0:bb[0],bb[1]],kidx)
   endif else begin
   kk=1.
   idxopoint=0
   endelse


   nisla=5
   dd=findgen(nisla)/(nisla-1)*abs((ee-aa))+ee
   print,'dd',dd,ee,ee0
   print,pror1(bb(0)),vectth(eeidx)
   print,'opoint',pror1(bb[0]+idxopoint(0)),vectth(bb[1])


;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;   print,'ee',ee
   phd_graphopener,filename=image_file,xsize=15.,ysize=10.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange=[0.,1.]
    yrange=[0.,2.*!pi]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1,vectth,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,/fill,col=0,yticks=5,ytickv=yt0,ytickname=ytpi,title=greekchi+'!Dhel!N '+strtrim(string(it,'(f9.4)')),xtitle='r/a',ytitle=greektheta
    contour,reform(rchi[*,*]),pror1,vectth,nlev=20,col=0,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,lev=ee[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1,vectth,levels=dd,col=23,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,levels=kk,col=43,thick=1.,/overplot

   phd_graphcloser,image_file
   image_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel_xy.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N '+strtrim(string(it,'(f9.4)')),xtitle='x/a',ytitle='y/a',/iso
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=223,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot

   phd_graphcloser,image_file
   image_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel_xy_fill.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N '+strtrim(string(it,'(f9.4)')),xtitle='x/a',ytitle='y/a',/iso,/fill,position=[0.20,0.15,0.94,0.94]
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=223,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot

   phd_graphcloser,image_file
;   endfor
  endelse



 
  
end
