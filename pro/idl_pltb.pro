pro idl_pltb,path,itp,m,n,d3,tok


;restore,path+'dat/spectrum.sav';,pror1,pror2,itp,bzt,btt,brt,mf_brt,mf_btt,mf_bzt
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
!p.charsize=1.5
; path='./'
; itp=3999
; m=1
; n=-10
 sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/bcyl.sav'
 ex_file = file_test(path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/bcyl.sav')   
 print,ex_file
 if (ex_file ne 1) then begin
    print,'need to create the magnetic field file for itp=',itp
;    call_procedure,'idl_saveb',path,itp
    restore,path+'dat/imf_bprofiles.sav'
    print,'max itp= ',n_elements(br[*,0,0,0])
    if (itp eq -1) then itp=n_elements(br[*,0,0,0])-1
    mf_brt = reform(br(itp,*,*,*))
    mf_btt = reform(bt(itp,*,*,*))
    mf_bzt = reform(bz(itp,*,*,*))
;#; check directory existence
    dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
    print,'dir',dir
    call_procedure,'check_dir',dir 
    sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/bcyl.sav'
    ttt=tt(itp) 
    save,filename=sav_file,pror1,pror2,ttt,mf_brt,mf_btt,mf_bzt
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
  jmn3=jmn3/m
  print,'tok',m,n,jmn3
  print, 'ATTENZIONE!!! RICONTROLLARE BENE SPETTRO!'
 endelse

help,mf_brt


 filename=path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/plot_m'+strtrim(string(m,format='(i0)'))+'_n'+strtrim(string(n,format='(i0)'))+'_b.eps'
 phd_graphopener,filename=filename,xsize=17.,ysize=10.

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter="164B;"
 greektau='!9' + String(letter) + '!X'
 greekphiletter="150B;"
 greekphi='!9' + String(greekphiletter) + '!X'

 xr = [0.,1.]
 if (tok ne 1) then begin
  yrb = [0.9*min(mf_brt[*,jmn,0]),max(mf_brt[*,jmn,0])*1.1]
  ytb = [0.9*min(mf_btt[*,jmn,0]),max(mf_btt[*,jmn,0])*1.1]
  yzb = [0.9*min(mf_bzt[*,jmn,0]),max(mf_bzt[*,jmn,0])*1.1]
 endif else begin
  yrb = [0.9*min(mf_brt[*,jmn3,0]),max(mf_brt[*,jmn3,0])*1.1]
  ytb = [0.9*min(mf_btt[*,jmn3,0]),max(mf_btt[*,jmn3,0])*1.1]
  yzb = [0.9*min(mf_bzt[*,jmn3,0]),max(mf_bzt[*,jmn3,0])*1.1]
 endelse
 phrng=[0.,2.*!pi]


; !p.charsize=0.7
; m1=0.08 & m2 =0.34 & m3=0.40 & m4=0.66 & m5 = 0.72 & m6 = 0.98
; m7=0.08 & m8=0.45 & m9 = 0.55 & m10 = 0.92
 m1=0.08 & m2 =0.33 & m3=0.41 & m4=0.66 & m5 = 0.73 & m6 = 0.98
 m7=0.08 & m8=0.45 & m9 = 0.55 & m10 = 0.92
 xyouts,0.34,0.955,'time = '+strtrim(string(ttt,format='(i0)'))+' '+greektau+'!DA!N, ITP='+strtrim(string(itp,'(i0)')),/norm,chars=1.2

; !p.charsize=1.2
;#; plot modbr
 !p.multi=[5,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yrb,xtitle='r/a',ytitle='mod b!Dr!N',title='mod b!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m9,m2,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yrb,col=0,  ytitle ='mod b!Dr!N',ystyle=1
 axis,yaxis=1,xrange=yrb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror1,mf_brt[*,jmn,0],thick=3.
 endif else begin 
 oplot,pror1,mf_brt[*,jmn3,0],thick=3.
 endelse

;#; plot phbr
 !p.multi=[4,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph b!Dr!N',title='phase b!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m7,m2,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph b!Dr!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
  oplot,pror1[1:-2],mf_brt[1:-2,jmn,1],thick=3.
 endif else begin 
  oplot,pror1[1:-2],mf_brt[1:-2,jmn3,1],thick=3.
 endelse

;#; plot modbt
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=ytb,xtitle='r/a',ytitle='mod b!D'+greektheta+'!N',title='mod b!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m9,m4,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=ytb,col=0,  ytitle ='mod b!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=ytb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2,mf_btt[*,jmn,0],thick=3.
 endif else begin 
 oplot,pror2,mf_btt[*,jmn3,0],thick=3.
 endelse
;#; plot phbt
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph b!D'+greektheta+'!N',title='phase b!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m7,m4,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph b!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2[1:-1],mf_btt[1:-1,jmn,1],thick=3.
 endif else begin 
 oplot,pror2[1:-1],mf_btt[1:-1,jmn3,1],thick=3.
 endelse

;#; plot modbz
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yzb,xtitle='r/a',ytitle='mod b!Dz!N',title='mod b!Dz!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m5,m9,m6,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yzb,col=0,  ytitle ='mod b!Dz!N',ystyle=1
 axis,yaxis=1,xrange=yzb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2,mf_bzt[*,jmn,0],thick=3.
 endif else begin 
 oplot,pror2,mf_bzt[*,jmn3,0],thick=3.
 endelse

;#; plot phbz
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph b!Dz!N',title='phase b!Dz!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m5,m7,m6,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph b!Dz!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 if (tok ne 1) then begin
 oplot,pror2[1:-1],mf_bzt[1:-1,jmn,1],thick=3.
 endif else begin 
 oplot,pror2[1:-1],mf_bzt[1:-1,jmn3,1],thick=3.
 endelse

 phd_graphcloser,filename

 print,'time is = ',ttt
end
