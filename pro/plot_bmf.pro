pro plot_bmf, path,d3,tmin,tmax,rds,tok


print,'ciao!'
letter="150B;"
greeketa='!9' + String(letter) + '!X'
letter="156B;"
greeknu='!9' + String(letter) + '!X'
letter="161B;"
greektheta='!9' + String(letter) + '!X'
loadct,0,file='~/idl/colors.tbl'
input='imf_bprofiles.sav'
restore,path+'dat/'+strtrim(input)

;result = file_test(path+'dat/spectrum.sav')
;if (result ne 1) then begin
; call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
; print,'nanf[1]',nanf[1]
;endif else begin
; restore,path+'dat/spectrum.sav'
; nzcon=nz
; nzcon[0]=nz[0]+1
;endelse
result = file_test(path+'dat/spectrum.dat')
if (result eq 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
 print,'nanf[1]',nanf[1],nz[0],nzcon[0]
endif else begin
 restore,path+'dat/spectrum.sav'
 nzcon=nz
 nzcon[0]=nz[0]+1
 if ((nzcon[0] - nz[0]) ge 2) then begin
  print,'spectrum.sav error in plot_bmf.pro'
  stop
 endif
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
if (tok eq 0) then begin
jin=69 ;usually corresponding to (1,-12)
jfin=76 ;usually corresponding to (1,-5)
;jfin=90 ;usually corresponding to (1,+9)
endif else begin 
jin=11
jfin= 23
endelse
endif else begin
jin=1
jfin=5
endelse

letter="164B;"
greektau='!9' + String(letter) + '!X'
dissipation = read_dissipation(path)
print,'dissipation',dissipation.eta,dissipation.nu
;#; Marco, mode naming
nmds = jfin - jin + 1
strings = make_array(nmds,/string)
colors = make_array(nmds,/integer)
help,strings
for k=0,nmds-1 do begin
 aa = mnum(jin+k,nzcon,nanf,d3)
 print,'prova',jin+k,aa
 strings(k) = strtrim('('+strtrim(string(aa(0),'(i0)'))+','+strtrim(string(aa(1),'(i0)'))+')')
 colors(k)=15+k*10
endfor
print,nzcon

if (d3 eq 1) then begin
 if (tok eq 0) then begin
  colors=[0,180,250,30,115,75,220,50]
  if (n_elements(colors) ne nmds) then stop
 endif else begin
  colors=[0,180,250,30,115,75,220,50,170,240,40,135,98]
  if (n_elements(colors) ne nmds) then stop
 endelse
endif
if (d3 eq 0) then begin
 colors=[20,180,250,90,115]
 if (n_elements(colors) ne nmds) then stop
endif
for k=0,nmds-1 do begin
 print,'modi & colori',jin+k,strings(k),colors(k)
endfor

print,nmds
help,colors
rds = float(rds)/100.
raggior = min(where(pror1 gt rds))
raggiot = min(where(pror2 gt rds))
print,'raggio=',raggior,raggiot

!p.multi=[0,1,2,0,0]
meld = indgen(jfin-jin+1) * 0. + 1.
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(br[tmin:tmax,raggior,jin:jfin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(br[see_mint:tmax,raggior,jin:jfin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
print,yrange
xrange=[min(tt),max(tt)]
help,tt
print,'tmax=',tt(tmax)
xrange=[tt(tmin),tt(tmax)]
print,'prova_xrange',xrange
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jin+2,0] / bt[tmin:tmax,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br) ',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 print,'bb',q,jin,colors(q-jin)
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
;#; marco for external nonresonant mode
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,87,0] / bt[tmin:tmax,101,0,0] *2.,col=0,thick=th*2.,psym=0
endif else begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.1
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(br)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 col=colors(q-jin)
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,1],col=col,thick=th,psym=0
 xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
endfor
print,nmds
for q=0,nmds-1l do begin
 col=colors(q)
 print,'qbr',colors(q),strings(q);,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor
; xyouts,xbase,ybase,'(1,-5)',col=50,/norm
; xyouts,xbase,ybase-deltay,'(1,+6)',col=0,/norm
; xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
; xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
; xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
; xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
; xyouts,xbase+3.*deltax,ybase,'(1,-11)',col=0,/norm
endif else begin
for q=jin,jfin do begin
oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,1],col=colors(q-jin),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+3.*deltax,ybase+deltay,'h='+strtrim(string(nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(3),/norm
endfor
endelse
xyouts,xbase+7.5*deltax,ybase+deltay,greeketa+'='+strtrim(string(float(dissipation.eta),format='(e11.4)')),col=0,/norm
xyouts,xbase+7.5*deltax,ybase+0.*deltay,greeknu+'='+strtrim(string(float(dissipation.nu),format='(e11.4)')),col=0,/norm

phd_graphcloser,filename

;#; Marco, temporaneo per figura EPS18
;filename=path+'dat/tempr_bprof_EPS18.eps'
;phd_graphopener,filename=filename,xsize=12,ysize=8
;
;plot,tt[tmin:tmax],br[tmin:tmax,raggior,jin+2,0] / bt[tmin:tmax,101,0,0] * 2.,position=[0.18,0.15,m1,0.85],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='mod(b!Dr!N) ',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
;if (d3 eq 1) then begin
;for q=jin,jfin do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
;endfor
;endif else begin
;for q=jin,jfin do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
;endfor
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,1,1] /(2.*!pi)*0.0065,col=colors(1-jin),thick=th,psym=0,linestyle=3
;endelse
;oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,1],col=colors(q-jin),thick=th,psym=0
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
 yrange=[0.,1.1*max(reform(bt[tmin:tmax,raggior,jin:jfin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(bt[see_mint:tmax,raggior,jin:jfin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jin+2,0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bt) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
endif else begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bt)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,1],col=colors(q-jin),thick=th,psym=0
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
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,1],col=colors(q-jin),thick=th,psym=0
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
 yrange=[0.,1.1*max(reform(bz[tmin:tmax,raggior,jin:jfin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(bz[see_mint:tmax,raggior,jin:jfin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jin+2,0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bz) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
endif else begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bz)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,1],col=colors(q-jin),thick=th,psym=0
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
for q=jin,jfin do begin
oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,1],col=colors(q-jin),thick=th,psym=0
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
meld = indgen(jfin-jin+1) * 0. + 1.
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(br[tmin:tmax,raggior,jin:jfin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(br[see_mint:tmax,raggior,jin:jfin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
tttt=indgen(fix(tmax)-fix(tmin))+fix(tmin)
;plot,tt[tmin:tmax],br[tmin:tmax,raggior,jin+2,0] / bt[tmin:tmax,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
plot,tttt,br[tmin:tmax,raggior,jin+2,0] / bt[tmin:tmax,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xtitle='',ytitle='mod(br) %',yr=yrange,ystyle=1,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
 oplot,tttt,br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
endif else begin
for q=jin,jfin do begin
; oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q),thick=th,psym=0
 oplot,tttt,br[tmin:tmax,raggior,q,0] / bt[tmin:tmax,101,0,0] *2.,col=colors(q-jin),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,tttt,br[tmin:tmax,raggior,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='ITP',ytitle='ph(br)',yr=[0.,2.*!pi],ystyle=1,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 col=colors(q-jin)
 oplot,tttt,br[tmin:tmax,raggior,q,1],col=col,thick=th,psym=0
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
for q=jin,jfin do begin
oplot,tttt,br[tmin:tmax,raggior,q,1],col=colors(q-jin),thick=th,psym=0
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
 oplot,tt[tmin:tmax],br[*,-1,0,0],col=0,thick=th,psym=0,linestyle=1
 oplot,tt[tmin:tmax],bt[*,-1,0,0],col=30,thick=th,psym=0,linestyle=1
 oplot,tt[tmin:tmax],bz[*,-1,0,0],col=240,thick=th,psym=0,linestyle=1


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],br[*,raggiot,0,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(b)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
 oplot,tt[tmin:tmax],br[*,raggiot,0,1],col=0,thick=th,psym=0
 oplot,tt[tmin:tmax],bt[*,raggiot,0,1],col=30,thick=th,psym=0
 oplot,tt[tmin:tmax],bz[*,raggiot,0,1],col=240,thick=th,psym=0
 oplot,tt[tmin:tmax],br[*,-1,0,1],col=0,thick=th,psym=0,linestyle=1
 oplot,tt[tmin:tmax],bt[*,-1,0,1],col=30,thick=th,psym=0,linestyle=1
 oplot,tt[tmin:tmax],bz[*,-1,0,1],col=240,thick=th,psym=0,linestyle=1
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

if (d3 eq 1) then begin ;#!# Marco, questo if si chiude alla fine
;#;Marco, modi m=0
 if (tok eq 0) then begin
  j0in = 17
  j0fin = 25

  nmds = j0fin - j0in + 1
  meld = indgen(j0fin-j0in+1) * 0. + 1.
  strings = make_array(nmds,/string)
  colors = make_array(nmds,/integer)
  help,strings
  for k=0,nmds-1 do begin
   aa = mnum(j0in+k,nzcon,nanf,d3)
   print,'prova m=0',j0in+k,aa
  strings(k) = strtrim('('+strtrim(string(aa(0),'(i0)'))+','+strtrim(string(aa(1),'(i0)'))+')')
  endfor

  if (d3 eq 1) then begin
   colors=[0,180,100,30,140,75,220,50,250]
   if (n_elements(colors) ne nmds) then stop
  endif
  if (d3 eq 0) then begin
   colors=[20,180,250,90,115]
   if (n_elements(colors) ne nmds) then stop
  endif
  for k=0,nmds-1 do begin
   print,'modi m0 & colori',j0in+k,strings(k),colors(k)
  endfor
 endif else begin ;#; modi m=0 per tokamak
  j0in = 1
  j0fin = 5
  nmds = j0fin - j0in + 1
  meld = indgen(j0fin-j0in+1) * 0. + 1.
  strings = make_array(nmds,/string)
  colors = make_array(nmds,/integer)
  help,strings
  for k=0,nmds-1 do begin
   aa = mnum(j0in+k,nzcon,nanf,d3)
   print,'prova m=0',j0in+k,aa
   strings(k) = strtrim('('+strtrim(string(aa(0),'(i0)'))+','+strtrim(string(aa(1),'(i0)'))+')')
  endfor
  colors=[20,180,250,90,115]
  for k=0,nmds-1 do begin
   print,'modi m0 & colori',j0in+k,strings(k),colors(k)
  endfor
 endelse

filename=path+'dat/tempt_m0_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(bt[tmin:tmax,raggiot,j0in:j0fin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(bt[see_mint:tmax,raggiot,j0in:j0fin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]

plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,j0in,0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(b!D'+greektheta+'!N)/B!D'+greektheta+'!N(a)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,j0in,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bt)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
 endfor
 for q=nmds-1l,0l,-1l do begin
  col=colors(q)
;  print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
  xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
  xyouts,xbase+3*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
  xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
  xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
  xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
 endfor
endelse

phd_graphcloser,filename

filename=path+'dat/tempr_m0_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(br[tmin:tmax,raggior,j0in:j0fin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(br[see_mint:tmax,raggior,j0in:j0fin,0] / bt[see_mint:tmax,101,0,0] *2.)]
 yrange=[0.,1.1*max(reform(br[see_mint:tmax,raggior,j0in:j0fin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]

plot,tt[tmin:tmax],br[tmin:tmax,raggior,j0in,0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(b!Dr!N)/B!D'+greektheta+'!N(a)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],br[tmin:tmax,raggior,j0in,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bt)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
 endfor
 for q=nmds-1l,0l,-1l do begin
  col=colors(q)
;  print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
  xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
  xyouts,xbase+3*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
  xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
  xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
  xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
 endfor
endelse

phd_graphcloser,filename

filename=path+'dat/tempz_m0_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(bz[tmin:tmax,raggiot,j0in:j0fin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(reform(bz[see_mint:tmax,raggiot,j0in:j0fin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]

plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,j0in,0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(b!Dz!N)/B!D'+greektheta+'!N(a)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,j0in,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bz)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
 endfor
 for q=nmds-1l,0l,-1l do begin
  col=colors(q)
;  print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
  xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
  xyouts,xbase+3*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
  xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
  xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
  xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
 endfor
endelse

phd_graphcloser,filename

filename=path+'dat/tempr2_m0_bprof.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
if (tmin ne 0.d) then begin
 yrange=[0.,1.1*max(reform(br[tmin:tmax,raggior,j0in:j0fin,0])/ reform(bt[tmin:tmax,101,0,0]#meld) *2.)]
endif else begin
 yrange=[0.,1.1*max(br[see_mint:tmax,raggior,j0in:j0fin,0] / bt[see_mint:tmax,101,0,0] *2.)]
 yrange=[0.,1.1*max(reform(br[see_mint:tmax,raggior,j0in:j0fin,0])/ reform(bt[see_mint:tmax,101,0,0]#meld) *2.)]
endelse
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
tttt=indgen(fix(tmax)-fix(tmin))+fix(tmin)

plot,tttt,br[tmin:tmax,raggior,j0in,0]/ bt[tmin:tmax,101,0,0] *2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(b!Dr!N)/B!D'+greektheta+'!N(a)',yr=yrange,ystyle=1,psym=0
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tttt,br[tmin:tmax,raggior,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tttt,br[tmin:tmax,raggior,q,0]/ bt[tmin:tmax,101,0,0] *2.,col=colors(q-j0in),thick=th,psym=0
 endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tttt,br[tmin:tmax,raggior,j0in,1],position=[m0,m2,m1,m3],/nodata,xtitle='time ('+greektau+'!DA!N)',ytitle='ph(bt)',yr=[0.,2.*!pi],ystyle=1
if (d3 eq 1) then begin
 for q=j0in,j0fin do begin
  oplot,tttt,br[tmin:tmax,raggior,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
 endfor
 for q=nmds-1l,0l,-1l do begin
  col=colors(q)
;  print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
  xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
 endfor
endif else begin
 for q=j0in,j0fin do begin
  oplot,tttt,br[tmin:tmax,raggior,q,1],col=colors(q-j0in),thick=th,psym=0
  xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
  xyouts,xbase+3*deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
  xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
  xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
  xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
 endfor
endelse

phd_graphcloser,filename


filename=path+'dat/tempr_m0_bprof.eps'
infiles = [path+'dat/tempr_m0_bprof.eps', path+'dat/tempt_m0_bprof.eps', path+'dat/tempz_m0_bprof.eps', path+'dat/tempr2_m0_bprof.eps']
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/m0_bprof.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor
endif;#; end if d3=1 for m=0 modes

end
