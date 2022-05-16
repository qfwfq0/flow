pro flow_xy,path,d3,tok,hel,itp0,itp1,deltaitp, mymin, mymax, nth, nzz,zed

;#; ANTIFFT parameters
; nth = 102
; nzz = 128

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter="167B;"
 greekomega='!9' + String(letter) + '!X'
 letter2="160B;"
 greekpi='!9' + String(letter2) + '!X'
 letter="143B;"
 greekchi='!9' + String(letter) + '!X'
;#; )check chosen times
 specyl_file = path+'dat/imf_bprofiles.sav'
 restore_bfile=path+'dat/ibprofiles.sav'
 part0 = test_up_to_date_files(restore_bfile,specyl_file)
 if (part0 eq 1) then begin
   call_procedure,'isavemodeprofiles',path,'0',d3
 endif
 restore,specyl_file
 if (itp1 eq -1) then begin 
  itp1=n_elements(br[*,0,0,0])-1
  itp0=n_elements(br[*,0,0,0])-1
  deltaitp=1
 endif
 if (itp0 eq itp1) then begin
;  outfile = path+'dat/itp/'+strtrim(string(itp0,format='(i0)'))+'/chihel.sav'
  print,'need to create the directory for helical flux function file for itp=',itp0
  dir =  path+'dat/itp/'+strtrim(string(itp0,format='(i0)'))+'/'
;  print,'dir',dir
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
;  print,'dir',dir
  call_procedure,'check_dir',dir 
 endfor
   
 print,'tok=',tok
 help,rchi
   ;path = './'
 if (rchi eq !null) then begin
;  file=path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/chihel.sav'
  inteli=fix(hel)
  file = path+'dat/itp/'+strtrim(string(itp0,format='(i0)'))+'/chi'+strtrim(string(inteli,format='(i0)'))+'.sav'
  fexis=file_test(file)
  print,'fexis',fexis
  if (fexis ne 1) then begin
;!  call_procedure, 'time_chihel',path,d3,tok,hel,itp0,itp1,deltaitp
;#; Marco, 14 marzo 2019: modo giusto per calcolare la funzione di flusso elicoidale
  print,itp0,itp1,deltaitp
  call_procedure, 'idl_savechi_v3',path,d3,tok,hel,nth,nzz,itp0,itp1,deltaitp
  endif
  restore,file
  vectth=th
  if (rchi2 ne !null) then begin
   rchi = reform(rchi2)
  endif else begin
  endelse
 endif

 if (vr eq !null) then begin
  file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/antifft_v_102x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  fexis=file_test(file)
  if (fexis ne 1) then begin
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
  endif
  restore,file
 endif
 if (br eq !null) then begin
  file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/antifft_b_102x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  fexis=file_test(file)
  if (fexis ne 1) then begin
   if (tok eq 1) then begin
    str='b'
    nth = nth
    nzz = nzz
    call_procedure, 'antifft_reconst_tok', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp, mymin, mymax
   endif else begin
    str='b'
    nth = nth
    nzz = nzz
;    call_procedure, 'antifft_reconst_rfp', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp, mymin, mymax
    call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
   endelse
  endif
  restore,file
 endif
;#; create current density
 if (jr eq !null) then begin
  file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/antifft_j_102x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  fexis=file_test(file)
  if (fexis ne 1) then begin
   if (tok eq 1) then begin
    str='j'
    nth = nth
    nzz = nzz
    call_procedure, 'antifft_reconst_tok', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp, mymin, mymax
   endif else begin
    str='j'
    nth = nth
    nzz = nzz
;    call_procedure, 'antifft_reconst_rfp', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp, mymin, mymax
    if (d3 ne 0) then begin
     call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
    endif
   endelse
  endif
  if (d3 ne 0) then begin
   restore,file
  endif
 endif
;#; create vorticity
 if (omr eq !null) then begin
  file = path+'dat/itp/'+strtrim(string(itp0,'(i0)'))+'/antifft_om_'+strtrim(string(102,'(i0)'))+'x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  print,'vorticity file',file,nth,nzz
  fexis=file_test(file)
  print,'exists vorticity file',fexis
  if (fexis ne 1) then begin
   if (tok eq 1) then begin
    str='om'
    nth = nth
    nzz = nzz
    call_procedure, 'antifft_reconst_tok', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp, mymin, mymax
   endif else begin
    str='om'
    nth = nth
    nzz = nzz
;    call_procedure, 'antifft_reconst_rfp', path, str, d3, hel ,nth, nzz, itp0, itp1, deltaitp, mymin, mymax
    if (d3 ne 0) then begin
     call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, itp0, itp1, deltaitp, mymin, mymax
    endif
   endelse
  endif
  if (d3 ne 0) then begin
   restore,file
  endif
 endif
 
 rchi = reform(rchi)
 help,vectth
 
 ppitp= indgen(nsaved)*deltaitp + itp0 
!p.charsize=1.5
 if (tok eq 1) then begin ;#; parte TOK
;#; ciclo sui tempi
 for pq=0,n_elements(ppitp)-1 do begin
  dir =  path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/'
  call_procedure,'check_dir',dir 

  str='b'
  outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  restore,outfile
  str='v'
  outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  restore,outfile
  str='j'
  outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(nth,format='(i0)'))+'x'+strtrim(string(nzz,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  restore,outfile
 ;#; restore helical flux function
  savechi_file = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/chi'+strtrim(string(inteli,format='(i0)'))+'.sav'
  restore,savechi_file

   image_file = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/flowxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
   ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]
   yt1=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]*4.
   ytpi1=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

   !p.psym=-3
;;#; interpolo br e vr su mesh x2
;   help,vr
;   help,br
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
     v2 = make_array(nn,nn,/double,value=0.d)
     xx = findgen(nn)/(nn-1)*2.-1.
     yy = findgen(nn)/(nn-1)*2.-1.
     dth = vecth(5) - vecth(4)
     dr = pror2(5) - pror2(4)
     for i=0,nn-1 do begin
      for k=0,nn-1 do begin
       if( (xx(i)^2.+yy(k)^2.) le 1. ) then begin
;#; calcolo r,th corrispondenti ai miei punti in x,y
        rr = sqrt(xx(i)^2.+yy(k)^2.)
        thv = atan(yy(k),xx(i))
        if (thv lt 0.d) then thv = thv + 2.*!pi
;#; calcolo gli indici frazionari per interpolare
        loc_th = value_locate(vecth,thv)
        idx_th = loc_th + (thv-vecth(loc_th))/dth
        loc_r = value_locate(pror2,rr)
        idx_r = loc_r + (rr-pror2(loc_r))/dr
;        print,idx_th,idx_r
;#; interpolo 
        tvr = interpolate(dvr[*,*,zzz],idx_r,idx_th,cubic=-0.5)
        tvt = interpolate(vt[*,*,zzz],idx_r,idx_th,cubic=-0.5)
        vx(i,k) = tvr*cos(thv) - tvt*sin(thv)
        vy(i,k) = tvr*sin(thv) + tvt*cos(thv)
       endif else begin
        vx(i,k) = 0.d
       endelse
        v2(i,k) = sqrt(vx(i,k)^2.+vy(i,k)^2.)
      endfor
     endfor
;#; velocity maximum
      maxv2 = max(v2)
     
;#; find separatrix
   if (vectth eq !null) then begin
    ;#; marco, array theta per la funzione di flusso elicoidale
    help,th
	stop
    vectth = th
   endif
   q=0
   aa = min(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   cc = max(rchi[bb[0],*]) ;#; separatrix!
;   print,'sep_lev=',cc
;   print,'Opoint_lev=',aa;,min(cc[0],aa[0])

   deltachi = max(rchi)-min(rchi)

   nisla=4
   dd=findgen(nisla)/(nisla-1)*abs((cc-aa))*0.75+aa
;   print,'dd',dd
;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;   print,'ee',ee
    print, 'image_file=',image_file
    help,tt
    print, 'AAA',ppitp(pq),tt(ppitp(pq))
    phd_graphopener,filename=image_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    print,'bbb',strtrim(string(tt[ppitp(pq)],'(f10.2)'))

    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot

    xyouts,0.6,0.85,/norm,'ITP'+strtrim(string(tt[ppitp(pq)],'(f10.2)')),chars=0.7
    xyouts,0.2,0.85,/norm,'max |v|'+strtrim(string(maxv2,'(e9.2)')),chars=0.7
    xyouts,0.6,0.9,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.7
    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
   phd_graphcloser,image_file


   loadct,13,/silent
   jr_file = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/jrxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jr_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(jr[*,*,0]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!Dr!N ITP'+strtrim(string(tt[ppitp(pq)],'(f10.2)')),xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(jr[*,*,0]),pror2#cos(vecth),pror2#sin(vecth),nlev=200,/fill,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
   phd_graphcloser,jr_file
   jt_file = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/jtxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jt_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!D'+greektheta+'!N ITP'+strtrim(string(tt[ppitp(pq)],'(f10.2)')),xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=200,/fill,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
   phd_graphcloser,jt_file
   jz_file = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/jzxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jz_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!Dz!N ITP'+strtrim(string(tt[ppitp(pq)],'(f10.2)')),xtitle='x/a',ytitle='y/a',/nodata,/iso
    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=200,/fill,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=184,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
   phd_graphcloser,jz_file

 infiles = make_array(3,/string)
 infiles = [jr_file,jt_file,jz_file]
 print, infiles
 infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/currentmyminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 ;print,infiles(i)
endfor

 endfor;#; fine del ciclo sui tempi: ATTENZIONE FARE ANCHE PER PARTE rfp
endif else begin ;#; end TOK part, begin RFP part, tok rfp
 
 for pq=0,n_elements(ppitp)-1 do begin
  print,'*******'
  print,'plotting time '+strtrim(string(ppitp(pq),'(i0)'))
  print,'*******'
  str='b'
  outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(102,format='(i0)'))+'x'+strtrim(string(nth,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  deltaa=1
  if (file_test(outfile) ne 1) then begin
   call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, ppitp(pq), ppitp(pq), deltaa, mymin, mymax
   restore,outfile
  endif else begin
   restore,outfile
  endelse
  str='v'
  outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(102,format='(i0)'))+'x'+strtrim(string(nth,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  if (file_test(outfile) ne 1) then begin
   call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, ppitp(pq), ppitp(pq), deltaa, mymin, mymax
   restore,outfile
  endif else begin
   restore,outfile
  endelse
  str='j'
  outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(102,format='(i0)'))+'x'+strtrim(string(nth,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  if (d3 ne 0) then begin
  if (file_test(outfile) ne 1) then begin
   call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, ppitp(pq), ppitp(pq), deltaa, mymin, mymax
   restore,outfile
  endif else begin
   restore,outfile
  endelse
  endif
  str='om'
  outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/antifft_'+str+'_'+strtrim(string(102,format='(i0)'))+'x'+strtrim(string(nth,format='(i0)'))+'myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.sav'
  if (d3 ne 0) then begin
  if (file_test(outfile) ne 1) then begin
   call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, ppitp(pq), ppitp(pq), deltaa, mymin, mymax
   restore,outfile
  endif else begin
   restore,outfile
  endelse
  endif
 ;#; restore helical flux function
  savechi_file = path+'dat/itp/'+strtrim(string(ppitp(pq),format='(i0)'))+'/chi'+strtrim(string(inteli,format='(i0)'))+'.sav'
  if (file_test(savechi_file) ne 1) then begin
   call_procedure, 'idl_savechi_v3',path,d3,tok,hel,nth,nzz,ppitp(pq),ppitp(pq),deltaa
  endif else begin
   restore,savechi_file
  endelse

  letter="161B;"
  greektheta='!9' + String(letter) + '!X'
  letter2="160B;"
  greekpi='!9' + String(letter2) + '!X'
  letter="143B;"
  greekchi='!9' + String(letter) + '!X'
  yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
  ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

  !p.psym=-3
;;; interpolo br e vr su mesh x2
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

;# scelgo il z e il theta
   zzz = zed
   ttt = 0
;# scrivo nuova mesh xx,yy
   nn = 51
   nn2 = 81
   nn3 = 31
   xx = make_array(nn,/double,value=0.d)
   yy = make_array(nn,/double,value=0.d)
   vx = make_array(nn,nn,/double,value=0.d)
   vy = make_array(nn,nn,/double,value=0.d)
   vx2 = make_array(nn3,nn2,/double,value=0.d)
   vzz = make_array(nn3,nn2,/double,value=0.d)
   xx = findgen(nn)/(nn-1)*2.-1.
   yy = findgen(nn)/(nn-1)*2.-1.
   xx2 = findgen(nn3)/(nn3-1)*2.-1.
   zz = findgen(nn2)/(nn2-1)*2.*!pi*4.
   dth = vecth(5) - vecth(4)
   dzz = vecz(5) - vecz(4)
   dr = pror2(5) - pror2(4)
   for i=0,nn-1 do begin
    for k=0,nn-1 do begin
     if( (xx(i)^2.+yy(k)^2.) le 1. ) then begin
;# alcolo r,th corrispondenti ai miei punti in x,y
      rr = sqrt(xx(i)^2.+yy(k)^2.)
      th = atan(yy(k),xx(i))
      if (th lt 0.d) then th = th + 2.*!pi
;# alcolo gli indici frazionari per interpolare
      loc_th = value_locate(vecth,th)
      idx_th = loc_th + (th-vecth(loc_th))/dth
      loc_r = value_locate(pror2,rr)
      idx_r = loc_r + (rr-pror2(loc_r))/dr
;      print,idx_th,idx_r
;# nterpolo 
      tvr = interpolate(dvr[*,*,zzz],idx_r,idx_th,cubic=-0.5)
      tvt = interpolate(vt[*,*,zzz],idx_r,idx_th,cubic=-0.5)
      vx(i,k) = tvr*cos(th) - tvt*sin(th)
      vy(i,k) = tvr*sin(th) + tvt*cos(th)
     endif else begin
      vx(i,k) = 0.d
     endelse
    endfor
   endfor

   for i=0,nn3-1 do begin
    for k=0,nn2-1 do begin
;# alcolo r,th corrispondenti ai miei punti in x,y
      rr = sqrt(xx2(i)^2.+yy(ttt)^2.)
      zzz = zz(k)
;# alcolo gli indici frazionari per interpolare
      loc_z = value_locate(vecz,zzz)
      idx_z = loc_z + (zzz-vecz(loc_z))/dzz
      loc_r = value_locate(pror2,rr)
      idx_r = loc_r + (rr-pror2(loc_r))/dr
;      print,idx_th,idx_r
;# nterpolo 
      caccc = reform(dvr[*,ttt,*])
      tvr2 = interpolate(caccc,idx_r,idx_z,cubic=-0.5,/grid)
      caccc = reform(vz[*,ttt,*])
      tvz = interpolate(caccc,idx_r,idx_z,cubic=-0.5)
      vx2(i,k) = tvr2
      vzz(i,k) = tvz
    endfor
   endfor
   
;#; find separatrix
;   aa = max(rchi[*,*],idx) ;#; O-point
;   bb = array_indices(rchi[*,*],idx)
;   ee0 = min(rchi[bb[0],*],eeidx)
;   ee = max(rchi[*,eeidx])
;   print,'sep_lev=',ee
;   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
;   deltachi = max(rchi)-min(rchi)

;#; find separatrix
   if (vectth eq !null) then begin
    vectth = vecth
   endif


   rchi = rchi[*,*,zed]
   aa = max(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   ee0 = min(rchi[bb[0],*],eeidx)
   ee = max(rchi[*,eeidx])
;   print,'bb',bb
;#; altra tecnica per trovare la separatrice guardando il minimo della derivata e poi cercando il massimo della chi a sinistra del minimo della derivata
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

;#; z section

   nisla=5
   dd=findgen(nisla)/(nisla-1)*abs((ee-aa))+ee
;   print,'dd',dd,ee,ee0
;   print,pror1(bb(0)),vectth(eeidx)
;   print,'opoint',pror1(bb[0]+idxopoint(0)),vectth(bb[1])
;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;   print,'ee',ee


    loadct,0,file='~/idl/colors.tbl'
   image_filexy = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/flowxy_vz_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_filexy,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/nodata,/iso
    nlev = 100
    colors = findgen(nlev)/(nlev-1)*252.
    if (n_elements(vecth) ne n_elements(reform(vz[0,*,zed]))) then begin
     vecth = findgen(n_elements(reform(vz[0,*,zed]))) / (n_elements(reform(vz[0,*,zed]))-1)*2.*!pi
    endif
    contour,reform(vz[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=nlev,c_colors=colors,title='v!Dz!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/
    nlev = 20
;    levs_vz = findgen(nlev)/(nlev-1) * 0.015 - 0.0075
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s0030p1000_0vms/s0030p1000_0vms_qinfty_fullspectrum/dat/itp
;    nlev = 20
;    levs_vz = findgen(nlev)/(nlev-1) * 0.0025 - 0.00125
; specifico per /ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/crash;
;    levs_vz = findgen(nlev)/(nlev-1) * 0.006 - 0.003
     levs_vz = findgen(nlev)/(nlev-1) * (max(vz[*,*,zed])-min(vz[*,*,zed])) + min(vz[*,*,zed])
    if (levs_vz eq !null) then begin
     contour,vz[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),nlev=nlev,/fill,/overplot
    endif else begin
     contour,vz[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_vz,/fill,/overplot
    endelse

;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=18,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee0[0],col=18,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot
;;#; marco to plot Opoint
;    plots,pror1[bb[0]+idxopoint[0]]*cos(vectth(bb[1])),pror1[bb[0]+idxopoint[0]]*sin(vectth(bb[1])),psym=symcat(16),col=0,syms=1.

    minv = min(vz[*,*,zed])
    maxv = max(vz[*,*,zed])

    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: v!Dz!N',chars=0.9
    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
    xyouts,0.1,0.84,/norm, 'max velocity '+strtrim(string(max(vz[*,*,zed]),'(f6.4)')),chars=0.9
    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
;    colorbar,/vertical,color=0,c_colors=colors,min=minv,max=maxv,position=[0.87, 0.11, 0.90, 0.90],/right,title='',charsize=charsize,format='(f7.3)',divisions=10
    if (levs_vz ne !null) then begin
    colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_vz),max=max(levs_vz),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
    endif
   phd_graphcloser,image_filexy

   image_filexyt = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/flowxy_vt_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_filexyt,xsize=12.,ysize=12.
    xrange=[-1.,1.]
;    xrange=[-0.3,0.3]

;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/
    nlev = 20
;    levs_vt = findgen(nlev)/(nlev-1) * 0.01 - 0.005
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s0030p1000_0vms/s0030p1000_0vms_qinfty_fullspectrum/dat/itp
;    nlev = 20
;    levs_vt = findgen(nlev)/(nlev-1) * 0.002 - 0.001
; specifico per /ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/crash;
;    levs_vt = findgen(nlev)/(nlev-1) * 0.006 - 0.003
; specifico per /ricercatori/ft/specyl/vivenzi/flow/archive/rfp/3d/eta_pol/prova6;
;    levs_vt = findgen(nlev)/(nlev-1) * 0.015 - 0.01
     levs_vt = findgen(nlev)/(nlev-1) * (max(vt[*,*,zed])-min(vt[*,*,zed])) + min(vt[*,*,zed])

    yrange = xrange
    contour,reform(vt[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,title='v!D'+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
    if (levs_vt eq !null) then begin
     contour,vt[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),nlev=nlev,/fill,/overplot
    endif else begin
     contour,vt[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_vt,/fill,/overplot
    endelse
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=18,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee0[0],col=18,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot

    minv = min(vt[*,*,zed])
    maxv = max(vt[*,*,zed])
    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: v!D'+greektheta+'!N',chars=0.9
    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
    xyouts,0.1,0.84,/norm, 'max velocity '+strtrim(string(max(vt[*,*,zed]),'(f6.4)')),chars=0.9
    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    if (levs_vt ne !null) then begin
    colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_vt),max=max(levs_vt),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
    endif
   phd_graphcloser,image_filexyt

   image_filexyr = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/flowxy_vr_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_filexyr,xsize=12.,ysize=12.
    xrange=[-1.,1.]
;    xrange=[-0.3,0.3]
    yrange = xrange
    contour,reform(vr[*,*,zed]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='v!Dr!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
    contour,vr[*,*,zed],pror1#cos(vecth),pror1#sin(vecth),nlev=20,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=18,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee0[0],col=18,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot

    minv = min(vr[*,*,zed])
    maxv = max(vr[*,*,zed])
    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: v!Dr!N',chars=0.9
    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
    xyouts,0.1,0.84,/norm, 'max velocity '+strtrim(string(max(vr[*,*,zed]),'(f6.4)')),chars=0.9
    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
   phd_graphcloser,image_filexyr

;#; velocity in the r-z plane
greekphiletter="146B;"
greekphi='!9' + String(greekphiletter) + '!X'
   image_filexz = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/flowxz_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_filexz,xsize=22.,ysize=22.
    yrange=[0.,2*!pi*4.]
    xrange = [-1.,1.]
    plot,[0],xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,xtitle='x/a',ytitle=greekphi,/nodata,position=[0.15,0.15,0.95,0.95]
    ;velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    velovect,vzz,vx2,xx2,zz,/overplot,col=20,min_value=1.e-5,length=1.5
   phd_graphcloser,image_filexz

infiles = [image_filexy,image_filexyt,image_filexyr,image_filexz]
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/flow'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
; spawn, "rm -f "+infiles(i)
; print,infiles(i)
endfor

 if (d3 ne 0) then begin
   image_filexy = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/bz_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_filexy,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    nlev = 100
    colors = findgen(nlev)/(nlev-1)*252.
    contour,reform(bz[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=nlev,c_colors=colors,title='B!Dz!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
; specifico per /ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/crash;
    levs_bz = findgen(nlev)/(nlev-1) * 1.8 - 0.9
    contour,bz[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_bz,/fill,/overplot
    minb = min(bz[*,*,zed])
    maxb = max(bz[*,*,zed])

    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: B!Dz!N',chars=0.9
;    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
;    xyouts,0.1,0.84,/norm, 'max velocity '+strtrim(string(max(vz[*,*,zed]),'(f6.4)')),chars=0.9
;    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_bz),max=max(levs_bz),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
   phd_graphcloser,image_filexy

   image_filexyt = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/bt_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_filexyt,xsize=12.,ysize=12.
    xrange=[-1.,1.]

; specifico per /ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/crash;
    levs_bt = findgen(nlev)/(nlev-1) * 1. - 0.5

    yrange = xrange
    contour,reform(bt[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,title='B!D'+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
    contour,bt[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_bt,/fill,/overplot

    minv = min(bt[*,*,zed])
    maxv = max(bt[*,*,zed])
    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: B!D'+greektheta+'!N',chars=0.9
;    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
;    xyouts,0.1,0.84,/norm, 'max velocity '+strtrim(string(max(vt[*,*,zed]),'(f6.4)')),chars=0.9
;    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_bt),max=max(levs_bt),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
   phd_graphcloser,image_filexyt

   image_filexyr = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/br_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=image_filexyr,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(br[*,*,zed]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='B!Dr!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
    contour,br[*,*,zed],pror1#cos(vecth),pror1#sin(vecth),nlev=20,/fill,/overplot

    minv = min(br[*,*,zed])
    maxv = max(br[*,*,zed])
    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: B!Dr!N',chars=0.9
;    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
    xyouts,0.1,0.84,/norm, 'max br '+strtrim(string(max(br[*,*,zed]),'(f6.4)')),chars=0.9
    xyouts,0.1,0.80,/norm, 'min br '+strtrim(string(min(br[*,*,zed]),'(f7.4)')),chars=0.9
;    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
   phd_graphcloser,image_filexyr

   loadct,13,/silent
   jr_file = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/jrxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jr_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(jr[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!Dr!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
    contour,reform(jr[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),nlev=200,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot
    minv = min(jr[*,*,zed])
    maxv = max(jr[*,*,zed])
    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: J!Dr!N',chars=0.9
    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
    xyouts,0.1,0.84,/norm, 'max current '+strtrim(string(max(jr[*,*,zed]),'(f6.4)')),chars=0.9
    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
    
   phd_graphcloser,jr_file
   jt_file = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/jtxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jt_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!D'+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
    contour,reform(jt[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=200,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot
    minv = min(jt[*,*,zed])
    maxv = max(jt[*,*,zed])
    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: J!D'+greektheta+'!N',chars=0.9
    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
    xyouts,0.1,0.84,/norm, 'max current '+strtrim(string(max(jt[*,*,zed]),'(f6.4)')),chars=0.9
    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
   phd_graphcloser,jt_file

   jz_file = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/jzxy_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
   phd_graphopener,filename=jz_file,xsize=12.,ysize=12.
    xrange=[-1.,1.]
    yrange = xrange
    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='J!Dz!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
    contour,reform(jz[*,*,0]),pror1#cos(vecth),pror1#sin(vecth),nlev=200,/fill,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=100,thick=2.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=25,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),nlev=20,col=0,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=184,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=254,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=254,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot
    minv = min(jz[*,*,zed])
    maxv = max(jz[*,*,zed])
    xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
    xyouts,0.1,0.97,/norm, 'filled contour: J!Dz!N',chars=0.9
    xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
    xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
    xyouts,0.1,0.84,/norm, 'max current '+strtrim(string(max(jz[*,*,zed]),'(f6.4)')),chars=0.9
    velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
    colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
   phd_graphcloser,jz_file


infiles = [jr_file,jt_file,jz_file]
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/currentmyminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
; print,infiles(i)
endfor
endif

;#; Marco, 26/04/2022: if d3=0 I don't compute vorticity
if (d3 ne 0) then begin
;#; plot vorticity
image_vortxyt = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/vort_vt_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
phd_graphopener,filename=image_vortxyt,xsize=12.,ysize=12.
xrange=[-1.,1.]
yrange = xrange
if (n_elements(vecth) ne n_elements(reform(omt[0,*,zed]))) then begin
 vecth = findgen(n_elements(reform(omt[0,*,zed]))) / (n_elements(reform(omt[0,*,zed]))-1)*2.*!pi
endif
contour,reform(omt[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title=greekomega+'!D'+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]

;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/
nlev = 20
levs_omt = findgen(nlev)/(nlev-1) * 0.2 - 0.1
levs_omt = findgen(nlev)/(nlev-1) * 0.09 - 0.045

levs_omt = findgen(nlev)/(nlev-1) * (max(omt[*,*,zed])-min(omt[*,*,zed])) + min(omt[*,*,zed])
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s0030p1000_0vms/s0030p1000_0vms_qinfty_fullspectrum/dat/itp
;nlev = 20
;levs_omt = findgen(nlev)/(nlev-1) * 0.06 - 0.03
contour,omt[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_omt,/fill,/overplot
;contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=250,thick=2.,/overplot
;contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=250,thick=1.,/overplot
levs_chi=reverse(-findgen(20)/19*0.2)
;contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=levs_chi,col=250,thick=2.,/overplot
minv = min(omt[*,*,zed])
maxv = max(omt[*,*,zed])
xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
xyouts,0.1,0.97,/norm, 'filled contour:'+greekomega+'!D'+greektheta+'!N',chars=0.9
xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
xyouts,0.1,0.85,/norm, 'max '+greekomega+' '+strtrim(string(max(omt[*,*,zed]),'(f6.4)')),chars=0.9
;velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_omt),max=max(levs_omt),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
phd_graphcloser,image_vortxyt

image_vortxyz = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/vort_vz_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
phd_graphopener,filename=image_vortxyz,xsize=12.,ysize=12.
xrange=[-1.,1.]
yrange = xrange
contour,reform(omz[*,*,zed]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,title=greekomega+'!Dz!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]

;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/
nlev = 20
levs_omz = findgen(nlev)/(nlev-1) * 0.2 - 0.1
;levs_omz = findgen(nlev)/(nlev-1) * 0.08 - 0.04
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s0030p1000_0vms/s0030p1000_0vms_qinfty_fullspectrum/dat/itp
;nlev = 20
;levs_omz = findgen(nlev)/(nlev-1) * 0.04 - 0.02
levs_omz = findgen(nlev)/(nlev-1) * (max(omz[*,*,zed])-min(omz[*,*,zed])) + min(omz[*,*,zed])

contour,omz[*,*,zed],pror1#cos(vecth),pror1#sin(vecth),lev=levs_omz,/fill,/overplot
;contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=250,thick=2.,/overplot
;contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=250,thick=1.,/overplot
;contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=levs_chi,col=250,thick=2.,/overplot
minv = min(omz[*,*,zed])
maxv = max(omz[*,*,zed])
xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
xyouts,0.1,0.97,/norm, 'filled contour:'+greekomega+'!Dz!N',chars=0.9
xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
xyouts,0.1,0.85,/norm, 'max '+greekomega+' '+strtrim(string(max(omz[*,*,zed]),'(f6.4)')),chars=0.9
;velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_omz),max=max(levs_omz),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
phd_graphcloser,image_vortxyz

image_vortxyr = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/vort_vr_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
phd_graphopener,filename=image_vortxyr,xsize=12.,ysize=12.
xrange=[-1.,1.]
yrange = xrange
contour,reform(omr[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title=greekomega+'!Dr!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]


levs_omr = findgen(nlev)/(nlev-1) * (max(omr[*,*,zed])-min(omr[*,*,zed])) + min(omr[*,*,zed])
;contour,omr[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),nlev=20,/fill,/overplot
contour,omr[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_omr,/fill,/overplot
contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=250,thick=2.,/overplot
contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=250,thick=1.,/overplot
contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=levs_chi,col=250,thick=2.,/overplot
minv = min(omr[*,*,zed])
maxv = max(omr[*,*,zed])
xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
xyouts,0.1,0.97,/norm, 'filled contour:'+greekomega+'!Dr!N',chars=0.9
xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
xyouts,0.1,0.85,/norm, 'max '+greekomega+' '+strtrim(string(max(omr[*,*,zed]),'(f6.4)')),chars=0.9
;velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
;colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_omr),max=max(levs_omr),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
phd_graphcloser,image_vortxyr

infiles = [image_vortxyr,image_vortxyt,image_vortxyz]
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/vort_v_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
; spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor
endif ;#; Marco, end if (d3 ne 0) for vorticity

;#; MARCO. qui faccio due plot dello flowshear
;#; check if files exists
file_shear=path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/shear_v_'+strtrim(string(102,'(i0)'))+'x'+strtrim(string(nth,'(i0)'))+'_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'
print,'shear file ',file_shear,' time: ',ppitp(pq)
if (file_test(file_shear) ne 1) then begin
 itp0 = ppitp(pq)
 itp1 = itp0
 delta = 1
 call_procedure,'flow_shear',path,d3,itp0,itp1,delta,nth,nzz,mymin,mymax
 restore,file_shear
endif else begin
 restore,file_shear
endelse

image_shearxyt = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/shear_vt_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
phd_graphopener,filename=image_shearxyt,xsize=12.,ysize=12.
xrange=[-1.,1.]
yrange = xrange
contour,reform(sh_t[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,title='shear v!D'+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/
nlev = 20
;levs_sht = findgen(nlev)/(nlev-1) * 0.1 - 0.05
;levs_sht = findgen(nlev)/(nlev-1) * 0.07 - 0.035
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s0030p1000_0vms/s0030p1000_0vms_qinfty_fullspectrum/dat/itp
;nlev = 20
;levs_sht = findgen(nlev)/(nlev-1) * 0.02 - 0.01
levs_sht = findgen(nlev)/(nlev-1) * (max(sh_t[*,*,zed])-min(sh_t[*,*,zed])) + min(sh_t[*,*,zed])

if (levs_sht eq !null) then begin
contour,sh_t[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),nlev=nlev,/fill,/overplot
endif else begin
contour,sh_t[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_sht,/fill,/overplot
endelse
minv = min(sh_t[*,*,zed])
maxv = max(sh_t[*,*,zed])
xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
xyouts,0.1,0.97,/norm, 'filled contour: shear v!D'+greektheta+'!N',chars=0.9
xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
xyouts,0.1,0.85,/norm, 'max shear '+strtrim(string(max(sh_t[*,*,zed]),'(f6.4)')),chars=0.9
;velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
if (levs_sht ne !null) then begin
colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_sht),max=max(levs_sht),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
endif
phd_graphcloser,image_shearxyt

image_shearxyz = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/shear_vz_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
phd_graphopener,filename=image_shearxyz,xsize=12.,ysize=12.
xrange=[-1.,1.]
yrange = xrange
contour,reform(sh_z[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,title='shear v!Dz!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s1000m10000/v2_nmp7_mp3e-3_qinfty/
nlev = 20
;levs_shz = findgen(nlev)/(nlev-1) * 0.2 - 0.1
;#; specifico per /0/rfx/ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/3d/theta16/s0030p1000_0vms/s0030p1000_0vms_qinfty_fullspectrum/dat/itp
;nlev = 20
;levs_shz = findgen(nlev)/(nlev-1) * 0.1 - 0.05
levs_shz = findgen(nlev)/(nlev-1) * (max(sh_z[*,*,zed])-min(sh_z[*,*,zed])) + min(sh_z[*,*,zed])

if (levs_shz eq !null) then begin
contour,sh_z[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),nlev=nlev,/fill,/overplot
endif else begin
contour,sh_z[*,*,zed],pror2#cos(vecth),pror2#sin(vecth),lev=levs_shz,/fill,/overplot
endelse
minv = min(sh_z[*,*,zed])
maxv = max(sh_z[*,*,zed])
xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
xyouts,0.1,0.97,/norm, 'filled contour: shear v!Dz!N',chars=0.9
xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
xyouts,0.1,0.85,/norm, 'max shear '+strtrim(string(max(sh_z[*,*,zed]),'(f6.4)')),chars=0.9
velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
if (levs_shz ne !null) then begin
colorbar_marco,/vertical,color=250,c_colors=colors,min=min(levs_shz),max=max(levs_shz),position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
endif
phd_graphcloser,image_shearxyz

image_shearxyr = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/shear_vr_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.eps'
phd_graphopener,filename=image_shearxyr,xsize=12.,ysize=12.
xrange=[-1.,1.]
yrange = xrange
contour,reform(sh_r[*,*,zed]),pror1#cos(vecth),pror1#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='shear v!Dr!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
contour,sh_r[*,*,zed],pror1#cos(vecth),pror1#sin(vecth),nlev=20,/fill,/overplot
minv = min(sh_r[*,*,zed])
maxv = max(sh_r[*,*,zed])
xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
xyouts,0.1,0.97,/norm, 'filled contour: shear v!Dr!N',chars=0.9
xyouts,0.1,0.93,/norm, 'arrows: vxy',chars=0.9
xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(ppitp(pq),'(i0)')),chars=0.9
xyouts,0.1,0.85,/norm, 'max shear '+strtrim(string(max(sh_r[*,*,zed]),'(f6.4)')),chars=0.9
velovect,vx,vy,xx,yy,/overplot,col=20,min_value=1.e-5,length=1.5
colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)';yticks=8,ytickv=[min,max/10000.,max/1000.,max/500.,max/100.,max/10.,max/5.,max/2.,max*1.];,ytickv=[strtrim(string(min,'(f5.3)')),strtrim(string(max/2,'(f5.3)')),strtrim(string(max,'(f5.3)'))];,title='!20L!Dc!N !M'+string(164B)+'!20a!7'
phd_graphcloser,image_shearxyr

infiles = [image_shearxyr,image_shearxyt,image_shearxyz]
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/itp/'+strtrim(string(ppitp(pq),'(i0)'))+'/shear_v_myminmax'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
; spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor


endfor ;#; end of cycle on time
endelse
end
