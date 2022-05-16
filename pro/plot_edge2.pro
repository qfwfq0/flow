pro plot_edge2, path,d3,tmin,tmax


loadct,0,file='~/idl/colors.tbl'
input='imf_bprofiles.sav'
restore,path+'dat/'+strtrim(input)
print,'restored!'


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
print,'tmin, tmax ',tmin,' ', tmax

;restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
filename=path+'dat/edgebr2.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88


colors=[0,250,20,115,75,220,50]
th=3.

;#; modi da plottare
if (d3 eq 1) then begin
jin=70 ;usually corresponding to (1,-11)
jfin=76 ;usually corresponding to (1,-5)
endif else begin
jin=1
jfin=2
endelse
rds = 1.d
raggior = -1
raggiot = -1
print,'raggio=',raggior,raggiot


!p.multi=[0,1,2,0,0]
yrange=[0.,1.2*max(br[tmin:tmax,raggior,jin:jfin,0])*2.]
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jin+2,0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 print,q,colors(q-jin)
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0]*2.,col=colors(q-jin),thick=th,psym=0
endfor
endif else begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],br[tmin:tmax,raggior,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='tau_A',ytitle='ph(br)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,1],col=colors(q-jin),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase,ybase,'(1,-5)',col=50,/norm
xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
endfor
endif else begin
for q=jfin,jin,-1 do begin
 oplot,tt[tmin:tmax],br[tmin:tmax,raggior,q,1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename
filename=path+'dat/edgebt2.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
yrange=[0.,1.2*max(bt[tmin:tmax,raggiot,jin:jfin,0])*2.]
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jin+2,0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bt)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,0]*2.,col=colors(q-jin),thick=th,psym=0
endfor
endif else begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bt[tmin:tmax,raggiot,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='tau_A',ytitle='ph(bt)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,1],col=colors(q-jin),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase,ybase,'(1,-5)',col=50,/norm
xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
endfor
endif else begin
for q=jfin,jin,-1 do begin
 oplot,tt[tmin:tmax],bt[tmin:tmax,raggiot,q,1],col=colors(q),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(2),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endfor
endelse

phd_graphcloser,filename
filename=path+'dat/edgebz2.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

!p.multi=[0,1,2,0,0]
yrange=[0.,1.2*max(bz[tmin:tmax,raggiot,jin:jfin,0])*2.]
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jin+2,0]*2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(bz)',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,0]*2.,col=colors(q-jin),thick=th,psym=0
endfor
endif else begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,0]*2.,col=colors(q),thick=th,psym=0
endfor
endelse


!p.multi=[1,1,2,0,0]
plot,tt[tmin:tmax],bz[tmin:tmax,raggiot,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='tau_A',ytitle='ph(bz)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,1],col=colors(q-jin),thick=th,psym=0
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggiot)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase,ybase,'(1,-5)',col=50,/norm
xyouts,xbase,ybase-deltay,'(1,-6)',col=200,/norm
xyouts,xbase+deltax,ybase,'(1,-7)',col=75,/norm
xyouts,xbase+deltax,ybase-deltay,'(1,-8)',col=115,/norm
xyouts,xbase+2.*deltax,ybase,'(1,-9)',col=20,/norm
xyouts,xbase+2.*deltax,ybase-deltay,'(1,-10)',col=250,/norm
endfor
endif else begin
for q=jfin,jin,-1 do begin
 oplot,tt[tmin:tmax],bz[tmin:tmax,raggiot,q,1],col=colors(q),thick=th,psym=0
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
infiles = [path+'dat/edgebr2.eps', path+'dat/edgebt2.eps', path+'dat/edgebz2.eps']
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/edgeb_zoom.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor
end
