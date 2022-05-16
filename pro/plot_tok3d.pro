pro plot_tok3d, path,d3,tmin,tmax,rds


letter="150B;"
greeketa='!9' + String(letter) + '!X'
letter="156B;"
greeknu='!9' + String(letter) + '!X'
loadct,0,file='~/idl/colors.tbl'
input='imf_bprofiles.sav'
restore,path+'dat/'+strtrim(input)
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
 nzcon=nz
 nzcon[0]=nz[0]-1
endelse

psym=0

print,keyword_set(tmax)
print,tmax
if (keyword_set(tmin) eq 0) then begin
 tmin = 0
 print,'not defined tmin', tmin
endif
if (keyword_set(tmax) eq 0) then begin
 tmax = n_elements(tt) - 1
 print,'not defined tmax', tmax
endif
if (tmax eq 1) then begin
 tmax = n_elements(tt) - 1
endif
print,'tmin, tmax ',tmin,' ', tmax,'   d3 ',d3
print,'aaa',tt(tmax)

filename=path+'dat/tempr_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88

if (n_elements(br[*,0,0]) gt 200) then begin
 see_mint = 50
endif else begin
 see_mint = 0
endelse

th=3.

;#; modi da plottare
if (d3 eq 1) then begin
j11=11 ;usually corresponding to (1,-1)
j12=22 ;usually corresponding to (2,-1)
j22=21;usually corresponding to (2,-2)
j32=32;usually corresponding to (3,-2)
j33=33;usually corresponding to (3,-3)
jplot=[j11,j12,j22,j32,j33]
endif else begin
jin=1
jfin=3
endelse

letter="164B;"
greektau='!9' + String(letter) + '!X'
dissipation = read_dissipation(path)
print,'dissipation',dissipation.eta,dissipation.nu
;#; Marco, mode naming
nmds = n_elements(jplot)
strings = make_array(nmds,/string)
colors = make_array(nmds,/integer)
help,strings
for k=0,nmds-1 do begin
 aa = mnum(jplot[k],nzcon,nanf,d3)
 strings(k) = strtrim('('+strtrim(string(aa(0),'(i0)'))+','+strtrim(string(aa(1),'(i0)'))+')')
 print,'prova',jplot(k),aa,strings(k)
endfor
;colors=[20,60,250,30,115,75,220,50]
colors=[20,60,250,230,115,75,30,50]

print,nmds
help,colors
rds = float(rds)/100.
raggior = min(where(pror1 gt rds))
raggiot = min(where(pror2 gt rds))
raggior=10
raggiot=10
print,'raggio=',raggior,raggiot

!p.multi=[0,1,2,0,0]
meld = indgen(nmds) * 0. + 1.
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(br[tmin:tmax,raggior,jplot,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(br[see_mint:tmax,raggior,jplot,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
print,'yrange=',yrange
;yrange=[0.,0.15]
xrange=[min(tt),max(tt)]
help,tt
print,'tmax=',tt(tmax)
xrange=[tt(tmin),tt(tmax)]
print,'prova_xrange',xrange
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[0],0] / bt[tmin:tmax,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br) ',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 print,'bb',q,colors(q)
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.1
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(br)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 col=colors(q)
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],1],col=col,thick=th,psym=0
 xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=0,nmds-1 do begin
 col=colors(q)
 print,'q',colors(q),strings(q);,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
; xyouts,xbase,ybase,'(1,-5)',col=50,/norm
; xyouts,xbase,ybase-deltay,'(1,-6)',col=220,/norm
; xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
; xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
; xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
; xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
; xyouts,xbase+3.*deltax,ybase,'(1,-11)',col=0,/norm
endif else begin
for q=0,nmds-1 do begin
oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+3.*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(3),/norm
endfor
endelse
xyouts,xbase+5*deltax,ybase+deltay,greeketa+'='+strtrim(string(float(dissipation.eta),format='(e11.4)')),col=0,/norm
xyouts,xbase+5*deltax,ybase+0.*deltay,greeknu+'='+strtrim(string(float(dissipation.nu),format='(e11.4)')),col=0,/norm

phd_graphcloser,filename

;#; Marco, temporaneo per figura EPS18
;filename=path+'dat/tempr_bprof_EPS18.eps'
;phd_graphopener,filename=filename,xsize=12,ysize=8
;
;plot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[0],0] / bt[tmin:tmax,101,0,0] * 2.,position=[0.18,0.15,m1,0.85],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='mod(b!Dr!N) ',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
;if (d3 eq 1) then begin
;for q=0,nmds-1 do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
;endfor
;endif else begin
;for q=0,nmds-1 do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
;endfor
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,1,1] /(2.*!pi)*0.0065,col=colors(1-jin),thick=th,psym=0,linestyle=3
;endelse
;oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],1],col=colors(q),thick=th,psym=0
;ybase=ybase-0.02
;deltay=deltay*1.5
;xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
;xyouts,xbase+3.*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
;xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
;xyouts,xbase+1.5*deltax,ybase,'(2,2h)',col=colors(1),/norm
;xyouts,xbase+3.*deltax,ybase,'(3,3h)',col=colors(2),/norm
;;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(3),/norm
;;endfor
;
;phd_graphcloser,filename
;stop
;
filename=path+'dat/tempt_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(bt[tmin:tmax,raggior,jplot,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(bt[see_mint:tmax,raggior,jplot,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[0],0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bt) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bt)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
;xyouts,xbase,ybase,'(1,-5)',col=50,/norm
;xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
;xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
;xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
;xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
;xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
;endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+3*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename

filename=path+'dat/tempz_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(bz[tmin:tmax,raggior,jplot,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(bz[see_mint:tmax,raggior,jplot,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[0],0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bz) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bz)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
;xyouts,xbase,ybase,'(1,-5)',col=50,/norm
;xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
;xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
;xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
;xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
;xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
;endfor
endif else begin
for q=0,nmds-1 do begin
oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+3.*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename

filename=path+'dat/tempr2_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15
!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(br[tmin:tmax,raggior,jplot,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(br[see_mint:tmax,raggior,jplot,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
;plot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[0],0] / bt[tmin:tmax,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
plot,br[tmin:tmax,raggior,jplot[0],0] / bt[tmin:tmax,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xtitle='',ytitle='mod(br) %',yr=yrange,ystyle=1,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
 oplot,br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
 oplot,br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,br[tmin:tmax,raggior,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='ITP',ytitle='ph(br)',yr=[0.,2.*!pi],ystyle=1,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 col=colors(q)
 oplot,br[tmin:tmax,raggior,jplot[q],1],col=col,thick=th,psym=0
 xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
; xyouts,xbase,ybase,'(1,-5)',col=50,/norm
; xyouts,xbase,ybase-deltay,'(1,-6)',col=220,/norm
; xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
; xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
; xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
; xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
; xyouts,xbase+3.*deltax,ybase,'(1,-11)',col=0,/norm
endif else begin
for q=0,nmds-1 do begin
oplot,br[tmin:tmax,raggior,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+3.*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse


phd_graphcloser,filename

filename=path+'dat/temp_bprof00.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


!p.multi=[0,1,2,0,0]
maxr=max(br[tmin:tmax,raggior,0,0])
maxt=max(bt[tmin:tmax,raggiot,0,0])
maxz=max(bz[tmin:tmax,raggiot,0,0])
maxx=max([maxr,maxt,maxz])
yrange=[0.,maxx*1.2]
plot,tt[tmin:tmax],br[*,raggiot,0,0],position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(b)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
 oplot,tt[tmin:tmax],br[*,raggiot,0,0],col=0,thick=th,psym=0
 oplot,tt[tmin:tmax],bt[*,raggiot,0,0],col=30,thick=th,psym=0
 oplot,tt[tmin:tmax],bz[*,raggiot,0,0],col=240,thick=th,psym=0


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],br[*,raggiot,0,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(b)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
 oplot,tt[tmin:tmax],br[*,raggiot,0,1],col=0,thick=th,psym=0
 oplot,tt[tmin:tmax],bt[*,raggiot,0,1],col=30,thick=th,psym=0
 oplot,tt[tmin:tmax],bz[*,raggiot,0,1],col=240,thick=th,psym=0

xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/lxx,format='(f6.3)')),col=0,/norm
xyouts,xbase,ybase,'br00',col=0,/norm
xyouts,xbase+deltax,ybase,'bt00',col=30,/norm
xyouts,xbase+2.*deltax,ybase,'bz00',col=240,/norm



phd_graphcloser,filename

filename=path+'dat/tempr_bprof.eps'
infiles = [path+'dat/tempr_bprof.eps', path+'dat/tempt_bprof.eps', path+'dat/tempz_bprof.eps', path+'dat/tempr2_bprof.eps', path+'dat/temp_bprof00.eps']
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/imf_bprof.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor

;#; VELOCITY PART
input='imf_vprofiles.sav'
restore,path+'dat/'+strtrim(input)

filename=path+'dat/tempr_vprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15
!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(vr[tmin:tmax,raggior,jplot,0])) *2.]
endif else begin
 yrange=[0.,1.1*max(reform(vr[see_mint:tmax,raggior,jplot,0])) *2.]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],vr[tmin:tmax,raggior,jplot[0],0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(vr)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vr[tmin:tmax,raggior,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vr[tmin:tmax,raggior,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03

lxx = n_elements(vr[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],vr[tmin:tmax,raggior,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(vr)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 col=colors(q)
 oplot,tt[tmin:tmax],vr[tmin:tmax,raggior,jplot[q],1],col=col,thick=th,psym=0
 xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
; xyouts,xbase,ybase,'(1,-5)',col=50,/norm
; xyouts,xbase,ybase-deltay,'(1,-6)',col=220,/norm
; xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
; xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
; xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
; xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
; xyouts,xbase+3.*deltax,ybase,'(1,-11)',col=0,/norm
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vr[tmin:tmax,raggior,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename
filename=path+'dat/tempt_vprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(vt[tmin:tmax,raggior,jplot,0])) *2.]
endif else begin
 yrange=[0.,1.1*max(reform(vt[see_mint:tmax,raggior,jplot,0])) *2.]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],vt[tmin:tmax,raggiot,jplot[0],0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(vt)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vt[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vt[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],vt[tmin:tmax,raggiot,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(vt)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vt[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
;xyouts,xbase,ybase,'(1,-5)',col=50,/norm
;xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
;xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
;xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
;xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
;xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
;endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vt[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+3*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename

filename=path+'dat/tempz_vprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(vz[tmin:tmax,raggior,jplot,0])) *2.]
endif else begin
 yrange=[0.,1.1*max(reform(vz[see_mint:tmax,raggior,jplot,0])) *2.]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],vz[tmin:tmax,raggiot,jplot[0],0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(vz)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vz[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vz[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],vz[tmin:tmax,raggiot,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(vz)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vz[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
;xyouts,xbase,ybase,'(1,-5)',col=50,/norm
;xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
;xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
;xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
;xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
;xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
;endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],vz[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename

filename=path+'dat/tempr2_vprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15
!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(vr[tmin:tmax,raggior,jplot,0])) *2.]
endif else begin
 yrange=[0.,1.1*max(reform(vr[see_mint:tmax,raggior,jplot,0])) *2.]
endelse
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
;plot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[0],0] / bt[tmin:tmax,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
plot,vr[tmin:tmax,raggior,jplot[0],0] * 2.,position=[m0,m4,m1,m5],/nodata,xtitle='ITP',ytitle='mod(vr) %',yr=yrange,ystyle=1,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
 oplot,vr[tmin:tmax,raggior,jplot[q],0] *2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
 oplot,vr[tmin:tmax,raggior,jplot[q],0]  *2.,col=colors(q),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03


!p.multi=[1,1,2,0,0]
plot,vr[tmin:tmax,raggior,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='ITP',ytitle='ph(vr)',yr=[0.,2.*!pi],ystyle=1,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 col=colors(q)
 oplot,vr[tmin:tmax,raggior,jplot[q],1],col=col,thick=th,psym=0
 xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
; xyouts,xbase,ybase,'(1,-5)',col=50,/norm
; xyouts,xbase,ybase-deltay,'(1,-6)',col=220,/norm
; xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
; xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
; xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
; xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
; xyouts,xbase+3.*deltax,ybase,'(1,-11)',col=0,/norm
endif else begin
for q=0,nmds-1 do begin
 oplot,vr[tmin:tmax,raggior,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse



phd_graphcloser,filename

filename=path+'dat/temp_vzprof00.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


!p.multi=[0,1,2,0,0]
if (tmin ne 0) then begin
maxr=max(vr[tmin:tmax,raggior,0,0])
maxt=max(vt[tmin:tmax,raggiot,0,0])
maxz=max(vz[tmin:tmax,raggiot,0,0])
maxr2=max(vr[tmin:tmax,1,0,0])
maxt2=max(vt[tmin:tmax,1,0,0])
maxz2=max(vz[tmin:tmax,*,0,0])
endif else begin
min_tt=0
maxr=max(vr[min_tt:tmax,raggior,0,0])
maxt=max(vt[min_tt:tmax,raggiot,0,0])
maxz=max(vz[min_tt:tmax,raggiot,0,0])
maxr2=max(vr[min_tt:tmax,*,0,0])
maxt3=max(vt[min_tt:tmax,10,0,0])
maxz3=max(vz[min_tt:tmax,10,0,0])
maxz2=max([maxz,maxz3])
maxt2=max([maxt,maxt3])
endelse
maxx=max([maxr,maxt,maxz,maxt2,maxz2])
yrange=[-maxz2*1.2,maxz2*1.2]
print,'CHECKrange, ',maxz2
if (maxz2 le 0.d) then begin
yrange=[-1.e-1,1.e-1]
endif
plot,tt[tmin:tmax],vr[*,raggiot,0,0],position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(vz)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
; oplot,tt[tmin:tmax],vr[*,raggiot,0,0],col=0,thick=th,psym=0
; oplot,tt[tmin:tmax],vt[*,raggiot,0,0],col=30,thick=th,psym=0
 oplot,tt[tmin:tmax],vz[*,raggiot,0,0]*cos(vz[*,raggiot,0,1]),col=240,thick=th,psym=0
; oplot,tt[tmin:tmax],vt[*,1,0,0],col=30,thick=th,linestyle=2,psym=0
 oplot,tt[tmin:tmax],vz[*,10,0,0]*cos(vz[*,10,0,1]),col=240,thick=th,linestyle=2,psym=0


!p.multi=[1,1,2,0,0]
yrange=[-maxx*1.1,maxx*1.1]
maxx=max([maxt,maxt3])

;yrange=[0.,maxx*1.2]
yrange=[-maxx*1.1,maxx*1.1]
if (maxx le 0.d) then begin
yrange=[-1.e-1,1.e-1]
endif
plot,tt[tmin:tmax],vr[*,raggiot,0,0]*cos(reform(vr[*,raggiot,0,1])),position=[m0,m2,m1,m3],/nodata,xtitle='time',ytitle='mod(vt)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
; oplot,tt[tmin:tmax],vr[*,raggiot,0,0]*cos(reform(vr[*,raggiot,0,1])),col=0,thick=th,psym=0
; oplot,tt[tmin:tmax],vt[*,raggiot,0,0]*cos(reform(vt[*,raggiot,0,1])),col=30,thick=th,psym=0
 oplot,tt[tmin:tmax],vt[*,raggiot,0,0]*cos(vt[*,raggiot,0,1]),col=20,thick=th,psym=0
; oplot,tt[tmin:tmax],vt[*,1,0,0]*cos(reform(vt[*,1,0,1])),col=30,thick=th,linestyle=2,psym=0
 oplot,tt[tmin:tmax],vt[*,10,0,0]*cos(vt[*,10,0,1]),col=20,thick=th,linestyle=2,psym=0
;plot,tt[tmin:tmax],vr[*,raggiot,0,1],position=[m0,m2,m1,m3],/nodata,xtitle='tau_A',ytitle='ph(v)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
; oplot,tt[tmin:tmax],vr[*,raggiot,0,1],col=0,thick=th,psym=0
; oplot,tt[tmin:tmax],vt[*,raggiot,0,1],col=30,thick=th,psym=0
; oplot,tt[tmin:tmax],vz[*,raggiot,0,1],col=240,thick=th,psym=0

xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03
xyouts,xbase,ybase+deltay,'SPC r/a='+strtrim(string(float(raggiot)/lxx,format='(f6.3)')),col=0,/norm
xyouts,xbase+1.5*deltax,ybase+deltay,'SPC r/a=0.1 (dashed)',col=0,/norm
xyouts,xbase,ybase,'vt00',col=20,/norm
;xyouts,xbase+deltax,ybase,'vt00',col=30,/norm
xyouts,xbase+2.*deltax,ybase,'vz00',col=240,/norm
xyouts,xbase+5*deltax,ybase+deltay,greeketa+'='+strtrim(string(float(dissipation.eta),format='(e11.4)')),col=0,/norm
xyouts,xbase+5*deltax,ybase+0.*deltay,greeknu+'='+strtrim(string(float(dissipation.nu),format='(e11.4)')),col=0,/norm

phd_graphcloser,filename


filename=path+'dat/temp_vrprof00.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

;if (tmin ne 0) then begin
;maxt=max(vt[tmin:tmax,raggiot,0,0])
;maxt2=max(vt[tmin:tmax,1,0,0])
;endif else begin
;maxt=max(vt[min_tt:tmax,raggiot,0,0])
;maxt2=max(vt[min_tt:tmax,1,0,0])
;endelse
;maxx=max([maxt,maxt2])
;yrange=[-maxx*1.1,maxx*1.1]

maxx=max([maxr,maxr2])
yrange=[0.,maxx*1.2]
yrange=[-maxx*1.1,maxx*1.1]
if (maxx le 0.d) then begin
yrange=[-1.e-1,1.e-1]
endif
plot,tt[tmin:tmax],vr[*,raggiot,0,0]*cos(reform(vr[*,raggiot,0,1])),position=[m0,m2,m1,m5],/nodata,xtitle='time',ytitle='mod(vr)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
 oplot,tt[tmin:tmax],vr[*,raggiot,0,0]*cos(reform(vt[*,raggiot,0,1])),col=30,thick=th,psym=0
 oplot,tt[tmin:tmax],vr[*,1,0,0]*cos(reform(vt[*,1,0,1])),col=30,thick=th,linestyle=2,psym=0
;plot,tt[tmin:tmax],vr[*,raggiot,0,1],position=[m0,m2,m1,m3],/nodata,xtitle='tau_A',ytitle='ph(v)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
; oplot,tt[tmin:tmax],vr[*,raggiot,0,1],col=0,thick=th,psym=0
; oplot,tt[tmin:tmax],vt[*,raggiot,0,1],col=30,thick=th,psym=0
; oplot,tt[tmin:tmax],vz[*,raggiot,0,1],col=240,thick=th,psym=0

xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03
xyouts,xbase,ybase+deltay,'SPC r/a='+strtrim(string(float(raggiot)/lxx,format='(f6.3)')),col=0,/norm
xyouts,xbase+1.5*deltax,ybase+deltay,'SPC r/a=0. (dashed)',col=0,/norm
xyouts,xbase+deltax,ybase,'vr00',col=30,/norm

phd_graphcloser,filename

filename=path+'dat/tempr_vprof.eps'
;infiles = [path+'dat/tempr_vprof.eps', path+'dat/tempt_vprof.eps', path+'dat/tempz_vprof.eps', path+'dat/tempr2_vprof.eps', path+'dat/temp_vzprof00.eps', path+'dat/temp_vtprof00.eps']
infiles = [path+'dat/temp_vzprof00.eps', path+'dat/tempr_vprof.eps', path+'dat/tempt_vprof.eps', path+'dat/tempz_vprof.eps', path+'dat/tempr2_vprof.eps']
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/imf_vprof.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor


;#; EDGE part


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88



rds = 1.d
raggior = -1
raggiot = -1

filename=path+'dat/edgebr.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15
!p.multi=[0,1,2,0,0]
yrange=[0.,1.2*max(br[tmin:tmax,raggior,jplot,0])*2.]
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[0],0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 print,q,colors(q)
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(br)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase,ybase,'(1,-5)',col=50,/norm
xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename
filename=path+'dat/edgebt.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
yrange=[0.,1.2*max(bt[tmin:tmax,raggiot,jplot,0])*2.]
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[0],0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bt)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bt)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase,ybase,'(1,-5)',col=50,/norm
xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename
filename=path+'dat/edgebz.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
yrange=[0.,1.2*max(bz[tmin:tmax,raggiot,jplot,0])*2.]
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[0],0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bz)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[0],1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bz)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase,ybase,'(1,-5)',col=50,/norm
xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
endfor
endif else begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jplot[q],1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse


phd_graphcloser,filename

filename=path+'dat/edgeb.eps'
infiles = [path+'dat/edgebr.eps', path+'dat/edgebt.eps', path+'dat/edgebz.eps']
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/edgeb.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor

end
