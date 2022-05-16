pro idl_pltb,path,itp,m,n

; path='./'
; itp=3999
; m=1
; n=-10
 sav_file = path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/bcyl.sav'
 ex_file = file_test(path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/bcyl.sav')   
 if (ex_file ne 1) then begin
    print,'need to create the magnetic field file for itp=',itp(q)
    call_procedure,'idl_saveb',path,itp
 endif
 restore,sav_file;,pror1,pror2,itp,bzt,btt,brt,mf_brt,mf_btt,mf_bzt
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
 jmn = janz2(m,n,my,mm,nanf,nz)



 filename=path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/plot_m'+strtrim(string(m,format='(i0)'))+'_n'+strtrim(string(n,format='(i0)'))+'_b.eps'
 phd_graphopener,filename=filename,xsize=17.,ysize=10.

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter="164B;"
 greektau='!9' + String(letter) + '!X'
 greekphiletter="150B;"
 greekphi='!9' + String(greekphiletter) + '!X'

 xr = [0.,1.]
 yrb = [0.9*min(mf_brt[0,*,jmn]),max(mf_brt[0,*,jmn])*1.1]
 ytb = [0.9*min(mf_btt[0,*,jmn]),max(mf_btt[0,*,jmn])*1.1]
 yzb = [0.9*min(mf_bzt[0,*,jmn]),max(mf_bzt[0,*,jmn])*1.1]
 phrng=[0.,2.*!pi]

 !p.charsize=0.7
 m1=0.06 & m2 =0.3 & m3=0.35 & m4=0.59 & m5 = 0.64 & m6 = 0.88
 m7=0.08 & m8=0.45 & m9 = 0.55 & m10 = 0.92
 xyouts,0.37,0.96,'time = '+strtrim(string(itp,format='(i0)'))+' '+greektau+'!DA!N',/norm

 !p.charsize=1.
;#; plot modbr
 !p.multi=[5,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yrb,xtitle='r',ytitle='mod b!Dr!N',title='mod b!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m9,m2,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yrb,col=0,  ytitle ='mod b!Dr!N',ystyle=1
 axis,yaxis=1,xrange=yrb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror1,mf_brt[0,*,jmn],thick=3.

;#; plot phbr
 !p.multi=[4,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r',ytitle='ph b!Dr!N',title='phase b!Dr!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m1,m7,m2,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph b!Dr!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror1,mf_brt[1,*,jmn],thick=3.

;#; plot modbt
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=ytb,xtitle='r',ytitle='mod b!D'+greektheta+'!N',title='mod b!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m9,m4,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=ytb,col=0,  ytitle ='mod b!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=ytb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_btt[0,*,jmn],thick=3.

;#; plot phbt
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r',ytitle='ph b!D'+greektheta+'!N',title='phase b!D'+greektheta+'!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m3,m7,m4,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph b!D'+greektheta+'!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_btt[1,*,jmn],thick=3.

;#; plot modbz
 !p.multi=[3,3,2,0,0]
 plot,[0.],xrange=xr,yrange=yzb,xtitle='r',ytitle='mod b!Dz!N',title='mod b!Dz!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m5,m9,m6,m10],ystyle=5,yminor=1,/nodata,chars=0.,xstyle=5
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=yzb,col=0,  ytitle ='mod b!Dz!N',ystyle=1
 axis,yaxis=1,xrange=yzb,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_bzt[0,*,jmn],thick=3.

;#; plot phbz
 !p.multi=[2,3,2,0,0]
 plot,[0.],xrange=xr,yrange=phrng,xtitle='r',ytitle='ph b!Dz!N',title='phase b!Dz!N m='+strtrim(string(m,format='(i0)'))+'n='+strtrim(string(n,format='(i0)')),position=[m5,m7,m6,m8],ystyle=5,yminor=1,/nodata,xstyle=5,chars=0.
 axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
 axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
 axis,yaxis=0,xrange=phrng,col=0,  ytitle ='ph b!Dz!N',ystyle=1
 axis,yaxis=1,xrange=phrng,col=0,  ytitle ='',ychars=0.01,ystyle=1
 oplot,pror2,mf_bzt[1,*,jmn],thick=3.


 phd_graphcloser,filename

end
