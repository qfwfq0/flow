pro plot_enbtok3d, path,d3,tmin,tmax

letter="164B;"
greektau='!9' + String(letter) + '!X'
letter="150B;"
greeketa='!9' + String(letter) + '!X'
letter="156B;"
greeknu='!9' + String(letter) + '!X'
!p.psym=0
loadct,0,file='~/idl/colors.tbl'
input=path+'dat/b_energy.sav'
input2=path+'dat/imf_bprofiles.sav'
ex_file = file_test(input)
print, 'ciao'
if (ex_file ne 1) then begin
 print,'need to create the magnetic field energy files'
 call_procedure,'idl_energy',path
endif

remake = test_up_to_date_files(input,input2)
if (remake ne 0) then begin
 print,'need to update the magnetic field energy files'
 call_procedure,'idl_energy',path
endif

restore,strtrim(input)
print,'restored!'

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
;tmin=0
;tmax=n_elements(tt)-1

dissipation = read_dissipation(path)
print,'dissipation',dissipation.eta,dissipation.nu

;restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
 nzcon=nz
 nz[0]=nz[0]-1
if (nzcon[0] eq nz[0]) then stop

endelse
filename=path+'dat/mag_en.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88

xbase=0.15
deltax=0.15
ybase=0.90
deltay=0.03


th=3.

;#; modi da plottare
if (d3 eq 1) then begin
j11=11 ;usually corresponding to (1,-1)
j12=22 ;usually corresponding to (2,-1)
j22=21;usually corresponding to (2,-2)
j32=32;usually corresponding to (3,-2)
j33=33;usually corresponding to (3,-3)
jplot=[j11,j12,j22,j32,j33]
nmds=n_elements(jplot)
endif else begin
jin=1
jfin=3
endelse

strings = make_array(nmds,/string)
if (d3 eq 1) then begin
colors=[20,60,250,230,115,75,30,50]
;if (n_elements(colors) ne nmds) then stop
endif
if (d3 eq 0) then begin
colors=[20,180,250,90,115]
;if (n_elements(colors) ne nmds) then stop
endif
for k=0,nmds-1 do begin
 aa = mnum(jplot[k],nzcon,nanf,d3)
 strings(k) = strtrim('('+strtrim(string(aa(0),'(i0)'))+','+strtrim(string(aa(1),'(i0)'))+')')
 print,'prova',jplot(k),aa,strings(k)
endfor
for k=0,nmds-1 do begin
 print,strings(k),colors(k)
endfor

help,b_comp
yrange=[0.,1.2*max(b_comp[tmin:tmax,jplot])]
print,'aaaaa',yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],b_comp[tmin:tmax,jplot[0]],/nodata,ytitle='W!DM!N',yr=yrange,ystyle=1,xr=xrange,xstyle=1
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],b_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
  print,'first',q,strings(q),colors(q)
endfor
for q=0,nmds-1 do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
  print,'second',strings(q),col
endfor
endif else begin
;#; 2d part
for q=0,nmds-1 do begin
 print,'first',q,colors(q)
 oplot,tt[tmin:tmax],b_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
endfor
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endelse
xyouts,0.7,ybase+deltay,greeketa+'='+strtrim(string(float(dissipation.eta),format='(e11.4)')),col=0,/norm
xyouts,0.7,ybase+0.*deltay,greeknu+'='+strtrim(string(float(dissipation.nu),format='(e11.4)')),col=0,/norm

phd_graphcloser,filename

;;#; modi da plottare
;if (d3 eq 1) then begin
;jin=70 ;usually corresponding to (1,-11)
;jfin=76 ;usually corresponding to (1,-5)
;endif else begin
;jin=0
;jfin=6
;endelse

filename=path+'dat/log_mag_en.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88


;colors=[0,250,20,115,75,220,50,0,0,0,0]
th=3.

help,b_comp
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
yrange=[1.e-6,1.e0]
yrange=[1.e-10,12*max(b_comp[tmin:tmax,jplot])]
plot,tt[tmin:tmax],b_comp[tmin:tmax,jplot[0]],/nodata,ytitle='W!DM!N',yr=yrange,ystyle=1,xr=xrange,xstyle=1,/ylog,xtitle='time ('+greektau+'!DA!N)',position=[0.11,0.11,0.94,0.94]
if (d3 eq 1) then begin
 oplot,tt[tmin:tmax],b_comp[tmin:tmax,0],col=0,thick=th
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],b_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor

endif else begin
;#; 2d part
 oplot,tt[tmin:tmax],b_comp[tmin:tmax,0],col=0,thick=th
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],b_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
endfor
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endelse
xyouts,0.7,ybase+deltay,greeketa+'='+strtrim(string(float(dissipation.eta),format='(e11.4)')),col=0,/norm
xyouts,0.7,ybase+0.*deltay,greeknu+'='+strtrim(string(float(dissipation.nu),format='(e11.4)')),col=0,/norm

phd_graphcloser,filename


input=path+'dat/v_energy.sav'
ex_file = file_test(input)
if (ex_file ne 1) then begin
 print,'need to create the kinetic field energy files'
 call_procedure,'idl_energy',path
endif
restore,strtrim(input)
print,'restored!'

;;#; modi da plottare
;if (d3 eq 1) then begin
;jin=70 ;usually corresponding to (1,-11)
;jfin=76 ;usually corresponding to (1,-5)
;endif else begin
;jin=1
;jfin=6
;endelse
filename=path+'dat/kin_en.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15

m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88


;colors=[0,250,20,115,75,220,50]
th=3.

help,v_comp
yrange=[0.,1.2*max(v_comp[tmin:tmax,jplot])]
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
plot,tt[tmin:tmax],v_comp[tmin:tmax,jplot[0]],/nodata,ytitle='W!DK!N',yr=yrange,ystyle=1,xr=xrange,xstyle=1
if (d3 eq 1) then begin
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],v_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
endfor
endif else begin
;#; oplot equilibrium part
oplot,tt[tmin:tmax],v_comp[tmin:tmax,0],col=0,thick=th*2.
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],v_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
endfor
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endelse

phd_graphcloser,filename

;#; modi da plottare
;if (d3 eq 1) then begin
;jin=70 ;usually corresponding to (1,-11)
;jfin=76 ;usually corresponding to (1,-5)
;endif else begin
;jin=0
;jfin=6
;endelse

filename=path+'dat/log_kin_en.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88


;colors=[0,250,20,115,75,220,50,0,0,0,0]
th=3.

help,v_comp
print,yrange
xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
yrange=[1.e-12,1.e-1]
plot,tt[tmin:tmax],v_comp[tmin:tmax,jplot[0]],/nodata,ytitle='W!DK!N',yr=yrange,ystyle=1,xr=xrange,xstyle=1,/ylog
if (d3 eq 1) then begin
 oplot,tt[tmin:tmax],v_comp[tmin:tmax,0],col=0,thick=th
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],v_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
endfor
endif else begin
oplot,tt[tmin:tmax],v_comp[tmin:tmax,0],col=0,thick=th*2.
for q=0,nmds-1 do begin
 oplot,tt[tmin:tmax],v_comp[tmin:tmax,jplot[q]],col=colors(q),thick=th
endfor
xbase = 0.55
deltax = 0.1
ybase = 0.3
deltay = 0.02
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors(0),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors(1),/norm
xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(2),/norm
;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm
endelse

phd_graphcloser,filename

filename=path+'dat/enb_spect.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88

th=3.

xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
yrange=[1.e-10,1.e-0]
;#; modes m=0
if (d3 ne 0) then begin
modes0 = [findgen(nz(0)-1)+nanf(0)]
modes1 = findgen(nz(1))+nanf(1)
modes2 = findgen(nz(2))+nanf(2)
endif else begin
modes0=0
modes1=1
modes2=2
endelse

xrange=[-30,10]
if (d3 ne 0) then begin
plot,modes0,b_comp[-1,0:nz(0)-1],/nodata,ytitle='W!DM!N',yr=yrange,ystyle=1,xr=xrange,xstyle=1,/ylog,xtitle='mode number'

oplot,modes0,b_comp[-1,1:nz(0)-1],psym=-symcat(16),syms=2.,linestyle=2,col=50
oplot,modes1,b_comp[-1,nz(0):nz(0)+nz(1)-1],psym=-symcat(16),syms=2.,linestyle=0,col=190
;oplot,modes2,b_comp[-1,nz(0)+nz(1):nz(0)+nz(1)+nz(2)-1],psym=-symcat(16),syms=1.,linestyle=3,col=240

endif else begin
endelse
;xbase = 0.5
;deltax = 0.1
;ybase = 0.3
;deltay = 0.02
;xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xbase=0.2
deltay=0.03
xyouts,xbase,ybase,'m=0',col=50,/norm
xyouts,xbase,ybase-deltay,'m=1',col=190,/norm
;xyouts,xbase,ybase-2.*deltay,'m=2',col=240,/norm
;xyouts,xbase+deltax,ybase,'(3,3h)',col=colors(3),/norm
;;xyouts,xbase+deltax,ybase-deltay,'(4,4h)',col=colors(4),/norm

phd_graphcloser,filename

filename=path+'dat/tot_en.eps'
phd_graphopener,filename=filename,xsize=20,ysize=15


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88


;colors=[0,250,20,115,75,220,50,0,0,0,0]
th=3.

xrange=[min(tt),max(tt)]
xrange=[tt(tmin),tt(tmax)]
yrange=[0.9*min(b_en+v_en),1.1*max(b_en+v_en)]/max(b_en+v_en)
plot,tt[tmin:tmax],b_en[tmin:tmax]+v_en[tmin:tmax],/nodata,ytitle='W!DK!N + W!DM!N',yr=yrange,ystyle=1,xr=xrange,xstyle=1,xtitle='time '+greektau+'!DA!N'
if (d3 eq 1) then begin
 oplot,tt[tmin:tmax],(b_en[tmin:tmax]+v_en[tmin:tmax])/max(b_en+v_en),col=0,thick=th*2.
endif else begin
 oplot,tt[tmin:tmax],(b_en[tmin:tmax]+v_en[tmin:tmax])/max(b_en+v_en),col=0,thick=th*2.
xbase = 0.5
deltax = 0.1
ybase = 0.3
deltay = 0.02
endelse

phd_graphcloser,filename



;#; create *pdf output
filename=path+'dat/mag_en.eps'
infiles = [path+'dat/log_mag_en.eps', path+'dat/mag_en.eps', path+'dat/log_kin_en.eps', path+'dat/kin_en.eps', path+'dat/tot_en.eps', path+'dat/enb_spect.eps']
infiles2 = strjoin(infiles, ' ')
outfile = path+'dat/en_prof.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor


end
