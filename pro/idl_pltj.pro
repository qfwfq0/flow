pro idl_pltj,path,itp,m,n,d3,tok


;restore,path+'dat/spectrum.sav';,pror1,pror2,itp,bzt,btt,brt,mf_jrt,mf_jtt,mf_jzt
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
 nzcon=nz
 nzcon[0]=nz[0]+1
endelse
if (d3 eq 0) then begin
 if (n ne nanf(m)) then begin
  print,'stop!, the mode is not in the 2D spectrum'
  exit
 endif
endif

!p.psym = 0
; path='./'
; itp=3999
; m=1
; n=-10
 ex_file = file_test(path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/jcyl.sav')   
 if (ex_file ne 1) then begin
    print,'need to create the current field file for itp=',itp
    file1=path+'dat/imf_jprofiles.sav'
    ex_file1 = file_test(file1)
    if (ex_file1 ne 1) then begin
      mf = 1
      call_procedure,'isavejprofiles',path,mf,d3
    endif
;    call_procedure,'idl_saveb',path,itp
    restore,path+'dat/imf_jprofiles.sav'
    if (itp eq -1) then itp=n_elements(mf_r[*,0,0,0])-1
    mf_jrt = reform(mf_r(itp,*,*,*))
    mf_jtt = reform(mf_t(itp,*,*,*))
    mf_jzt = reform(mf_z(itp,*,*,*))
    pror1 = findgen(n_elements(mf_t[0,*,0,0]))/(n_elements(mf_t[0,*,0,0]) - 1)
    pror2 = findgen(n_elements(mf_r[0,*,0,0]))/(n_elements(mf_t[0,*,0,0]) - 1) - 1.d/(2.d*n_elements(mf_t[0,*,0,0]))
;#; check directory existence
    dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
    sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/jcyl.sav'
    print,'dir',dir
    call_procedure,'check_dir',dir 
    ttt=t(itp) 
    save,filename=sav_file,pror1,pror2,ttt,mf_jrt,mf_jtt,mf_jzt
 endif
 sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/jcyl.sav'
 print,sav_file
 restore,sav_file;,pror1,pror2,itp,bzt,btt,brt,mf_jrt,mf_jtt,mf_jzt
 if (tok ne 1) then begin
  jmn = janz2(m,n,my,mm,nanf,nzcon)
  print,'check_RFP',m,n,jmn
 endif else begin
  jmn = janz3(m,n,my,mm,nanf,nz)
  print,'tok',m,n,jmn
 endelse

help,mf_jrt
print,jmn
print,'ttt= ',ttt


 filename=path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/plot_m'+strtrim(string(m,format='(i0)'))+'_n'+strtrim(string(n,format='(i0)'))+'_j.eps'
 phd_graphopener,filename=filename,xsize=17.,ysize=10.

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter="164B;"
 greektau='!9' + String(letter) + '!X'
 greekphiletter="150B;"
 greekphi='!9' + String(greekphiletter) + '!X'

 xr = [0.,1.]
 yrb = [0.9*min(mf_jrt[*,jmn,0]),max(mf_jrt[*,jmn,0])*1.1]
 ytb = [0.9*min(mf_jtt[*,jmn,0]),max(mf_jtt[*,jmn,0])*1.1]
 yzb = [0.9*min(mf_jzt[*,jmn,0]),max(mf_jzt[*,jmn,0])*1.1]
 phrng=[0.,2.*!pi]

 !p.charsize=0.7
 m1=0.06 & m2 =0.32 & m3=0.38 & m4=0.64 & m5 = 0.70 & m6 = 0.96
 m7=0.08 & m8=0.45 & m9 = 0.55 & m10 = 0.92
 xyouts,0.37,0.96,'time = '+strtrim(string(ttt,format='(i0)'))+' '+greektau+'!DA!N',/norm,chars=1.3

 !p.charsize=1.2
;#; plot modbr
 !p.multi=[5,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yrb,xtitle='r/a',ytitle='mod j!Dr!N',title='mod j!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m9,m2,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yrb,col=0,  ytitle ='mod j!Dr!N',ystyle=1
 axis,yaxis=1,xrange=yrb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror1,mf_jrt[*,jmn,0],thick=3.

;#; plot phbr
 !p.multi=[4,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph j!Dr!N',title='phase j!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m7,m2,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph j!Dr!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror1,mf_jrt[*,jmn,1],thick=3.

;#; plot modbt
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=ytb,xtitle='r/a',ytitle='mod j!D'+greektheta+'!N',title='mod j!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m9,m4,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=ytb,col=0,  ytitle ='mod j!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=ytb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_jtt[*,jmn,0],thick=3.

;#; plot phbt
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph j!D'+greektheta+'!N',title='phase j!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m7,m4,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph j!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_jtt[*,jmn,1],thick=3.

;#; plot modbz
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yzb,xtitle='r/a',ytitle='mod j!Dz!N',title='mod j!Dz!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m5,m9,m6,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yzb,col=0,  ytitle ='mod j!Dz!N',ystyle=1
 axis,yaxis=1,xrange=yzb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_jzt[*,jmn,0],thick=3.

;#; plot phbz
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r/a',ytitle='ph j!Dz!N',title='phase j!Dz!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m5,m7,m6,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph j!Dz!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_jzt[*,jmn,1],thick=3.


 phd_graphcloser,filename
 print,'time is = ',ttt

end
