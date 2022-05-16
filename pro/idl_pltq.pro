pro idl_pltq,path,itp,d3,tok


result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
; call_procedure, 'read_spectrum',path,my,mm,nanf,nz
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
;restore,path+'dat/spectrum.sav';,pror1,pror2,itp,bzt,btt,brt,mf_brt,mf_btt,mf_bzt
;if (d3 eq 0) then begin
; if (n ne nanf(m)) then begin
;  print,'stop!, the mode is not in the 2D spectrum'
;  exit
; endif
;endif

!p.psym = 0
; path='./'
; itp=3999
; m=1
; n=-10
qfile=path+'dat/qprof.sav'
 ex_file = file_test(path+'dat/qprof.sav')   
 if (ex_file ne 1) then begin
    print,'need to create the safety factor file'
    call_procedure,'qprof',path
    restore,qfile
    if (itp eq -1) then itp=n_elements(qprof[*,0])-1
;#; check directory existence
    dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
    print,'dir',dir
    call_procedure,'check_dir',dir 
    ttt=tt(itp) 
 endif

 aaa = test_up_to_date_files(qfile,path+'dat/imf_bprofiles.sav')
 if (aaa eq 1) then begin
  call_procedure,'qprof',path
  restore,qfile
 endif

sav_file=path+'dat/qprof.sav'
restore,sav_file;,pror1,pror2,itp,bzt,btt,brt,mf_brt,mf_btt,mf_bzt
if (itp eq -1) then itp=n_elements(qprof[*,0])-1
dir =  path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/'
print,'dir',dir
call_procedure,'check_dir',dir 
if (itp eq -1) then itp=n_elements(qprof[*,0])-1

ttt=tt(itp) 
print,'time=',ttt

if (tok ne 1) then begin
 filename=path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/qplot.eps'
 phd_graphopener,filename=filename,xsize=7.,ysize=5.

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter="164B;"
 greektau='!9' + String(letter) + '!X'
 greekphiletter="150B;"
 greekphi='!9' + String(greekphiletter) + '!X'

 xr = [0.,1.]
 if (min(qprof[itp,*]) gt 0.d) then begin
 yr = [0.9*min(qprof[itp,*]),1.1*max(qprof[itp,*])]
 endif else begin
 yr = [1.1*min(qprof[itp,*]),1.1*max(qprof[itp,*])]
 endelse
 
; yr=[-0.05,0.20]
 phrng=[0.,2.*!pi]


 !p.charsize=0.7
 m1=0.18 & m2 =0.34 & m3=0.40 & m4=0.66 & m5 = 0.72 & m6 = 0.96
 m7=0.2 & m8=0.45 & m9 = 0.55 & m10 = 0.91
; xyouts,0.37,0.9,'time = '+strtrim(string(ttt,format='(f9.0)'))+' '+greektau+'!DA!N',/norm,chars=0.8

print,'yr_rfp',yr
;print,m1,m7,m6,m10
;#; plot modbr
; plot,[0.],xrange=xr,yrange=yr,xtitle='r/a',ytitle='q(r/a)',title='safety factor',position=[m1,m7,m6,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
; plot,[0.],xrange=xr,yrange=yr,xtitle='r/a',ytitle='q(r/a)',title='safety factor',/nodata
; axis,xaxis=0,xrange=xr,col=0,  xtitle ='r/a'
; axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
; axis,yaxis=0,xrange=yr,col=0,  ytitle ='q(r/a)',ystyle=1
; axis,yaxis=1,xrange=yr,col=0,  ytitle ='',ychars=0.01,ystyle=1
 plot,pror2,qprof[itp,*],thick=3.,position=[m1,m7,m6,m10],xtitle='r/a',ytitle='q(r/a)',xstyle=1,ystyle=1,xrange=xr,yrange=yr,/nodata,title='axisymmetric q'
 oplot,pror2,qprof[itp,*],col=0,thick=3
 plots,xr[0],0.
 plots,xr[1],0.,/continue,linestyle=2,thick=2.
 loadct,12,/silent
 plots,xr[0],1./6.
 plots,xr[1],1./6.,/continue,linestyle=2,thick=2.,col=30
 xo=0.8 
 yo = 0.65
 dyo=0.054
 xyouts,xo,yo,'q=1/6',/norm,chars=0.8,col=30
 loadct,0,file='~/idl/colors.tbl'
 plots,xr[0],1./7.
 plots,xr[1],1./7.,/continue,linestyle=2,thick=2.,col=220
 plots,xr[0],1./8.
 plots,xr[1],1./8.,/continue,linestyle=2,thick=2.,col=90
 plots,xr[0],1./9.
 plots,xr[1],1./9.,/continue,linestyle=2,thick=2.,col=20
; plots,xrange[0],1./9.
; plots,xrange[1],1./9.,/continue,linestyle=2,thick=2.,col=0
 xyouts,xo,yo-1.*dyo,'q=1/7',/norm,chars=0.8,col=220
 xyouts,xo,yo-2.*dyo,'q=1/8',/norm,chars=0.8,col=90
 xyouts,xo,yo-3.*dyo,'q=1/9',/norm,chars=0.8,col=20
 phd_graphcloser,filename

endif else begin ;#; tokamak case

 filename=path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/qplot.eps'
 phd_graphopener,filename=filename,xsize=7.,ysize=5.

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter="164B;"
 greektau='!9' + String(letter) + '!X'
 greekphiletter="150B;"
 greekphi='!9' + String(greekphiletter) + '!X'

 xr = [0.,1.]
 yr = [0.9*min(qprof[itp,*]),1.1*max(qprof[itp,*])]
 phrng=[0.,2.*!pi]

 !p.charsize=0.7
 m1=0.16 & m2 =0.34 & m3=0.40 & m4=0.66 & m5 = 0.72 & m6 = 0.96
 m7=0.2 & m8=0.45 & m9 = 0.55 & m10 = 0.92

 plot,pror2,qprof[itp,*],thick=3.,position=[m1,m7,m6,m10],xtitle='r/a',ytitle='q(r/a)',xstyle=1,ystyle=1,xrange=xr,yrange=yr,/nodata
 oplot,pror2,qprof[itp,*],col=0,thick=3
 plots,xr[0],1.
 plots,xr[1],1.,/continue,linestyle=2,thick=2.,col=30
 xyouts,0.4,0.65,'q=1',/norm,chars=0.8,col=30
 plots,xr[0],2.
 plots,xr[1],2.,/continue,linestyle=2,thick=2.,col=220
 plots,xr[0],3.
 plots,xr[1],3.,/continue,linestyle=2,thick=2.,col=190
 xyouts,0.4,0.7,'q=2',/norm,chars=0.8,col=220
 xyouts,0.4,0.75,'q=3',/norm,chars=0.8,col=190
 h=1
 res = min(where(qprof[itp,*] gt 1./h))
 plots,pror2(res),yr[0]
 plots,pror2(res),yr[1],/continue,linestyle=1,thick=1.
 h=0.5
 res = min(where(qprof[itp,*] gt 1./h))
 plots,pror2(res),yr[0]
 plots,pror2(res),yr[1],/continue,linestyle=1,thick=1.
 phd_graphcloser,filename

;#; plotto q(0) in time
 filename=path+'dat/time_q0.eps'
 phd_graphopener,filename=filename,xsize=7.,ysize=5.

 xr = [0.,max(tt)]
 yr = [0.9*min(q02),1.1*max(q95)]
 phrng=[0.,2.*!pi]

 !p.charsize=0.7
 m1=0.16 & m2 =0.34 & m3=0.40 & m4=0.66 & m5 = 0.72 & m6 = 0.96
 m7=0.2 & m8=0.45 & m9 = 0.55 & m10 = 0.92

 plot,tt,q02,thick=3.,position=[m1,m7,m6,m10],xtitle='time'+greektau+'!DA!N',ytitle='q(0)',xstyle=1,ystyle=1,xrange=xr,yrange=yr,/nodata
 oplot,tt,q02,col=0,thick=3
 oplot,tt,q95,col=210,thick=3
 plots,xr[0],1.
 plots,xr[1],1.,/continue,linestyle=2,thick=2.,col=30
 xyouts,0.4,0.65,'q=1',/norm,chars=0.8,col=30
 plots,xr[0],2.
 plots,xr[1],2.,/continue,linestyle=2,thick=2.,col=220
 plots,xr[0],3.
 plots,xr[1],3.,/continue,linestyle=2,thick=2.,col=190
 xyouts,0.4,0.7,'q=2',/norm,chars=0.8,col=220
 xyouts,0.4,0.75,'q=3',/norm,chars=0.8,col=190
; h=1
; res = min(where(qprof[itp,*] gt 1./h))
; plots,pror2(res),yr[0]
; plots,pror2(res),yr[1],/continue,linestyle=1,thick=1.
; h=0.5
; res = min(where(qprof[itp,*] gt 1./h))
; plots,pror2(res),yr[0]
; plots,pror2(res),yr[1],/continue,linestyle=1,thick=1.
 phd_graphcloser,filename
endelse

 print,'time is = ',ttt
end
