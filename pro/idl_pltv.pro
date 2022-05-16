pro idl_pltv,path,itp,m,n,d3,tok


;restore,path+'dat/spectrum.sav';,pror1,pror2,itp,bzt,btt,brt,mf_brt,mf_btt,mf_bzt
print,path
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
 nzcon=nz
 nzcon[0]=nz[0]+1
endelse
m=fix(m)
n=fix(n)
if (d3 eq 0) then begin
 if (n ne nanf(m)) then begin
  print,'stop!, the mode is not in the 2D spectrum',n,m,nanf
  exit
 endif
endif

!p.psym = 0
; path='./'
; itp=3999
; m=1
; n=-10
 sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/vcyl.sav'
 ex_file = file_test(path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/vcyl.sav')   
 if (ex_file ne 1) then begin
    print,'need to create the velocity field file for itp=',itp
;    call_procedure,'idl_saveb',path,itp
    restore,path+'dat/imf_vprofiles.sav'
    print,'max itp= ',n_elements(vr[*,0,0,0])
    if (itp eq -1) then itp=n_elements(vr[*,0,0,0])-1
    mf_vrt = reform(vr(itp,*,*,*))
    mf_vtt = reform(vt(itp,*,*,*))
    mf_vzt = reform(vz(itp,*,*,*))
;#; check directory existence
    dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
    print,'dir',dir
    call_procedure,'check_dir',dir 
    ttt=tt(itp) 
    sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/vcyl.sav'
    save,filename=sav_file,pror1,pror2,ttt,mf_vrt,mf_vtt,mf_vzt
 endif
 restore,sav_file;,pror1,pror2,itp,bzt,btt,brt,mf_brt,mf_btt,mf_bzt
 if (itp eq -1) then itp=n_elements(br[*,0,0,0])-1
 dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
 call_procedure,'check_dir',dir 
 if (tok ne 1) then begin
  jmn = janz2(m,n,my,mm,nanf,nzcon)
  print,'check_RFP',m,n,jmn
 endif else begin
  jmn3 = janz3(m,n,my,mm,nanf,nz)
  print,'tok',m,n,jmn3
 endelse

help,mf_brt

print,'time=',ttt

 filename=path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/plot_m'+strtrim(string(m,format='(i0)'))+'_n'+strtrim(string(n,format='(i0)'))+'_v.eps'
 phd_graphopener,filename=filename,xsize=17.,ysize=10.

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter="164B;"
 greektau='!9' + String(letter) + '!X'
 greekphiletter="150B;"
 greekphi='!9' + String(greekphiletter) + '!X'

 xr = [0.,1.]
 if (tok ne 1) then begin
 yrb = [0.9*min(mf_vrt[*,jmn,0]),max(mf_vrt[*,jmn,0])*1.1]
 ytb = [0.9*min(mf_vtt[*,jmn,0]),max(mf_vtt[*,jmn,0])*1.1]
 yzb = [0.9*min(mf_vzt[*,jmn,0]),max(mf_vzt[*,jmn,0])*1.1]
 endif else begin
 yrb = [0.9*min(mf_vrt[*,jmn3,0]),max(mf_vrt[*,jmn3,0])*1.1]
 ytb = [0.9*min(mf_vtt[*,jmn3,0]),max(mf_vtt[*,jmn3,0])*1.1]
 yzb = [0.9*min(mf_vzt[*,jmn3,0]),max(mf_vzt[*,jmn3,0])*1.1]
 endelse
 phrng=[0.,2.*!pi]


 !p.charsize=0.7
 m1=0.08 & m2 =0.33 & m3=0.41 & m4=0.66 & m5 = 0.74 & m6 = 0.985
 m7=0.08 & m8=0.45 & m9 = 0.55 & m10 = 0.92
 xyouts,0.37,0.96,'time = '+strtrim(string(ttt,format='(f9.0)'))+' '+greektau+'!DA!N',/norm,chars=1.3

 !p.charsize=1.2
;#; plot modbr
 !p.multi=[5,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yrb,xtitle='r/a',ytitle='mod v!Dr!N',title='mod v!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m9,m2,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yrb,col=0,  ytitle ='mod v!Dr!N',ystyle=1
 axis,yaxis=1,xrange=yrb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror1,mf_vrt[*,jmn,0],thick=3.
 endif else begin
 oplot,pror1,mf_vrt[*,jmn3,0],thick=3.
 endelse

;#; plot phbr
 !p.multi=[4,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph v!Dr!N',title='phase v!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m7,m2,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph v!Dr!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror1[1:-2],mf_vrt[1:-2,jmn,1],thick=3.
 endif else begin
 oplot,pror1[1:-2],mf_vrt[1:-2,jmn3,1],thick=3.
 endelse

;#; plot modbt
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=ytb,xtitle='r/a',ytitle='mod v!D'+greektheta+'!N',title='mod v!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m9,m4,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=ytb,col=0,  ytitle ='mod v!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=ytb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2,mf_vtt[*,jmn,0],thick=3.
 endif else begin
 oplot,pror2,mf_vtt[*,jmn3,0],thick=3.
 endelse

;#; plot phbt
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph v!D'+greektheta+'!N',title='phase v!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m7,m4,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph v!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2[1:-2],mf_vtt[1:-2,jmn,1],thick=3.
 endif else begin
 oplot,pror2[1:-2],mf_vtt[1:-2,jmn3,1],thick=3.
 endelse

;#; plot modbz
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yzb,xtitle='r/a',ytitle='mod v!Dz!N',title='mod v!Dz!N m='+strtrim(string(m,format='(i0)'))+' n='+strtrim(string(n,format='(i0)')),position=[m5,m9,m6,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yzb,col=0,  ytitle ='mod v!Dz!N',ystyle=1
 axis,yaxis=1,xrange=yzb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2,mf_vzt[*,jmn,0],thick=3.
 endif else begin
 oplot,pror2,mf_vzt[*,jmn3,0],thick=3.
 endelse

;#; plot phbz
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph v!Dz!N',title='phase v!Dz!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m5,m7,m6,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph v!Dz!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2[1:-2],mf_vzt[1:-2,jmn,1],thick=3.
 endif else begin
 oplot,pror2[1:-2],mf_vzt[1:-2,jmn3,1],thick=3.
 endelse


 phd_graphcloser,filename

 print,'time is = ',ttt
end
