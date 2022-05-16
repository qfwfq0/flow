pro compare, path0, path1, d3,tmin,tmax

;path0 = '/ricercatori/ft/specyl/veranda/flow/archive/rfp/2d/theta16/momsour/h10/s300m30_momsour2e-3/'
str00 = strpos(path0,'rfp')
str0 = strmid(path0,str00)
l0 = strlen(str0)
st0 = strmid(str0,0,l0-1)
;path1 = '/ricercatori/ft/specyl/veranda/flow/archive/rfp/2d/theta16/momsour/h10/s300m30_momsour0/'
str11 = strpos(path1,'rfp')
str1 = strmid(path1,str11)
l1 = strlen(str1)
st1 = strmid(str1,0,l1-1)

remchar,st0,'/'
remchar,st1,'/'

loadct,0,file='~/idl/colors.tbl'
input='imf_bprofiles.sav'
print,path0
if (br0 eq !null) then begin
restore,path0+'dat/'+strtrim(input)
br0 = br
bt0 = bt
bz0 = bz
tt0 = tt
result = file_test(path0+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path0,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path0+'dat/spectrum.sav'
endelse
print,mm
print,nanf
print,nz
restore,path1+'dat/'+strtrim(input)
br1 = br
bt1 = bt
bz1 = bz
tt1 = tt
endif
input2='imf_vprofiles.sav'
print,path0
if (vr0 eq !null) then begin
restore,path0+'dat/'+strtrim(input2)
vr0 = vr
vt0 = vt
vz0 = vz
tt0 = tt
result = file_test(path1+'dat/spectrum.sav')
if (result ne 1) then begin
; call_procedure, 'read_spectrum',path1,my,mm,nanf,nz,nzcon
 call_procedure, 'read_spectrum',path1,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path1+'dat/spectrum.sav'
endelse
restore,path1+'dat/'+strtrim(input2)
vr1 = vr
vt1 = vt
vz1 = vz
tt1 = tt
endif
print,'restored!'

psym=0

print,keyword_set(tmax)
;print,tmax
;if (keyword_set(tmin) eq 0) then begin
 tmin = 0
 print,'not defined tmin', tmin
;endif
;if (keyword_set(tmax) eq 0) then begin
 tmax0 = n_elements(tt0) - 1
 tmax1 = n_elements(tt1) - 1
 tminmax = min([tt0(tmax0),tt1(tmax1)])
 idxmax0 = min(where(tt0 gt tminmax))
 idxmax1 = min(where(tt1 gt tminmax))
;endif



filename=path0+'magn_compare_'+st0+'-'+st1+'.eps'
print,filename
phd_graphopener,filename=filename,xsize=20,ysize=15


m0=0.12
m1=0.92
m2=0.12
m3=0.48
m4=0.52
m5=0.88


colors=[10,180,250,30,115,75,220,50]
colors0=[10,180,250,30,115,75,220,50]
colors=[180,250,30,115,75,220,50]
colors0=[180,250,30,115,75,220,50]
colors1=reverse(colors0)
colors1=reverse(colors0)
th=3.

;#; modi da plottare
if (d3 eq 1) then begin
jin=70 ;usually corresponding to (1,-11)
jfin=76 ;usually corresponding to (1,-5)
endif else begin
jin=1
jfin=2
endelse

;#; Marco, mode naming
nmds = jfin - jin + 1
strings = make_array(nmds,/string)
help,strings
for k=0,nmds-1 do begin
 aa = mnum(jin+k,nz,nanf,d3)
 strings(k) = strtrim('('+strtrim(string(aa(0),'(i0)'))+','+strtrim(string(aa(1),'(i0)'))+')')
endfor

rds = 0.70
raggior = min(where(pror1 gt rds))
raggiot = min(where(pror2 gt rds))
print,'raggio=',raggior,raggiot

!p.multi=[0,1,2,0,0]
meld = indgen(jfin-jin+1) * 0. + 1.
if (tmin ne 0.d) then begin
 bm0=  max(reform(br0[tmin:idxmax0,raggior,jin:jfin,0])/ reform(bt0[tmin:idxmax0,101,0,0]#meld) *2.)
 bm1=  max(reform(br1[tmin:idxmax1,raggior,jin:jfin,0])/ reform(bt1[tmin:idxmax1,101,0,0]#meld) *2.)
 bmm = max([bm0,bm1])
 yrange = [0,bmm]
endif else begin
 bm0 = max(reform(br0[50:idxmax0,raggior,jin:jfin,0])/ reform(bt0[50:idxmax0,101,0,0]#meld) *2.)
 bm1 = max(reform(br1[50:idxmax1,raggior,jin:jfin,0])/ reform(bt1[50:idxmax1,101,0,0]#meld) *2.)
 bmm = max([bm0,bm1])
 yrange = [0,bmm]
endelse
print,yrange
xrange=[min(tt),tminmax]
print,'prova_xrange',xrange
plot,tt0[tmin:idxmax0],br0[tmin:idxmax0,raggior,jin+2,0] / bt0[tmin:idxmax0,101,0,0] * 2.,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(br) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt0[tmin:idxmax0],br0[tmin:idxmax0,raggior,q,0] / bt0[tmin:idxmax0,101,0,0] *2.,col=colors0(q-jin),thick=th,psym=0
 oplot,tt1[tmin:idxmax1],br1[tmin:idxmax1,raggior,q,0] / bt1[tmin:idxmax1,101,0,0] *2.,col=colors0(q-jin),thick=th,psym=0,linestyle=2
endfor
endif else begin
for q=jin,jfin do begin
 oplot,tt0[tmin:idxmax0],br0[tmin:idxmax0,raggior,q,0] / bt0[tmin:idxmax0,101,0,0] *2.,col=colors0(q),thick=th,psym=0
 oplot,tt1[tmin:idxmax1],br1[tmin:idxmax1,raggior,q,0] / bt1[tmin:idxmax1,101,0,0] *2.,col=colors0(q),thick=th,psym=0,linestyle=2
endfor
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
plot,tt0[tmin:idxmax0],br0[tmin:idxmax0,raggior,jin+2,1],position=[m0,m2,m1,m3],/nodata,xtitle='tau_A',ytitle='ph(br)',yr=[0.,2.*!pi],ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt0[tmin:idxmax0],br0[tmin:idxmax0,raggior,q,1],col=colors0(q-jin),thick=th,psym=0
 oplot,tt1[tmin:idxmax1],br1[tmin:idxmax1,raggior,q,1],col=colors0(q-jin),thick=th,psym=0,linestyle=2
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
xyouts,xbase+0.5*deltax,ybase-15.6*deltay,st0,col=0,/norm
xyouts,xbase+0.5*deltax,ybase-14.6*deltay,st1,col=0,/norm
xyouts,xbase+0.1*deltax,ybase-14.6*deltay,'dashed ',col=250,/norm
endif else begin
for q=jfin,jin,-1 do begin
 oplot,tt0[tmin:idxmax0],br0[tmin:idxmax0,raggior,q,1],col=colors0(q),thick=th,psym=0
 oplot,tt1[tmin:idxmax1],br1[tmin:idxmax1,raggior,q,1],col=colors1(q),thick=th,psym=0,linestyle=2
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors0(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors0(2),/norm
xyouts,xbase+0.6*deltax,ybase,'(3,3h)',col=colors0(3),/norm
xyouts,xbase+0.6*deltax,ybase-deltay,'(4,4h)',col=colors0(4),/norm

endfor
jin=1
xyouts,xbase+1.1*deltax,ybase,st0,col=colors0(1),/norm
xyouts,xbase+1.1*deltax,ybase-deltay,st1,col=colors1(1),/norm
endelse

phd_graphcloser,filename

filename_v=path0+'vel_compare_'+st0+'-'+st1+'.eps'
print,filename_v
phd_graphopener,filename=filename_v,xsize=20,ysize=15

;#; Marco, mode naming
raggior = min(where(pror1 gt rds))
raggiot = min(where(pror2 gt rds))
print,'raggio=',raggior,raggiot

!p.multi=[0,1,2,0,0]
meld = indgen(jfin-jin+1) * 0. + 1.
if (tmin ne 0.d) then begin
 vm0=  max(reform(vz0[tmin:idxmax0,raggior,0,0]))
 vm1=  max(reform(vz1[tmin:idxmax1,raggior,0,0]))
 vmm = max([vm0,vm1])
 yrange = [0,vmm]
endif else begin
 vm0 = max(reform(vz0[50:idxmax0,raggior,0,0]))
 vm1 = max(reform(vz1[50:idxmax1,raggior,0,0]))
 vmm = max([vm0,vm1])
 yrange = [0,vmm]
endelse
print,yrange
xrange=[min(tt),tminmax]
print,'prova_xrange',xrange,yrange
plot,tt0[tmin:idxmax0],vz0[tmin:idxmax0,raggior,0,0] ,position=[m0,m4,m1,m5],/nodata,xchars=0.01,ytitle='mod(vz00) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
 oplot,tt0[tmin:idxmax0],vz0[tmin:idxmax0,raggior,0,0],col=0,thick=th,psym=0
 oplot,tt1[tmin:idxmax1],vz1[tmin:idxmax1,raggior,0,0],col=250,thick=th,linestyle=2,psym=0
endif else begin
 oplot,tt0[tmin:idxmax0],vz0[tmin:idxmax0,raggior,0,0] ,col=0,thick=th,psym=0
 oplot,tt1[tmin:idxmax1],vz1[tmin:idxmax1,raggior,0,0] ,col=250,thick=th,linestyle=2,psym=0
endelse


xbase=0.10
deltax=0.2
ybase=0.93
deltay=0.03

lxx = n_elements(br[0,*,0,0]) - 1

!p.multi=[1,1,2,0,0]
if (tmin ne 0.d) then begin
 vm0=  max(reform(vt0[tmin:idxmax0,raggior,0,0]))
 vm1=  max(reform(vt1[tmin:idxmax1,raggior,0,0]))
 vmm = max([vm0,vm1])
 yrange = [0,vmm]
endif else begin
 vm0 = max(reform(vt0[50:idxmax0,raggior,0,0]))
 vm1 = max(reform(vt1[50:idxmax1,raggior,0,0]))
 vmm = max([vm0,vm1])
 yrange = [0,vmm]
endelse
print,yrange
xrange=[min(tt),tminmax]
print,'prova_xrange',xrange,yrange
plot,tt0[tmin:idxmax0],vt0[tmin:idxmax0,raggior,0,0],position=[m0,m2,m1,m3],/nodata,ytitle='mod(vt00) %',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0
if (d3 eq 1) then begin
 oplot,tt0[tmin:idxmax0],vt0[tmin:idxmax0,raggior,0,0],col=0,thick=th,psym=0
 oplot,tt1[tmin:idxmax1],vt1[tmin:idxmax1,raggior,0,0],col=250,thick=th,linestyle=2,psym=0
endif else begin
 oplot,tt0[tmin:idxmax0],vt0[tmin:idxmax0,raggior,0,0] ,col=0,thick=th,psym=0
 oplot,tt1[tmin:idxmax1],vt1[tmin:idxmax1,raggior,0,0] ,col=250,thick=th,linestyle=2,psym=0
endelse

;xyouts,xbase+0.4*deltax,ybase,st0,col=1,/norm
;xyouts,xbase+0.4*deltax,ybase-deltay,st1,col=0,/norm
;xyouts,xbase+0.0*deltax,ybase-deltay,'dashed ',col=250,/norm
xyouts,xbase+0.5*deltax,ybase-15.6*deltay,st0,col=0,/norm
xyouts,xbase+0.5*deltax,ybase-14.6*deltay,st1,col=250,/norm
xyouts,xbase+0.1*deltax,ybase-14.6*deltay,'dashed ',col=250,/norm
raggiorr=pror2(raggior)
xyouts,xbase+0.1*deltax,ybase,'SPC r/a ='+strtrim(string(raggiorr,'(f6.3)')),col=0,/norm
phd_graphcloser,filename_v

input='b_energy.sav'
print,path0
if (b_comp0 eq !null) then begin
restore,path0+'dat/'+strtrim(input)
b_comp0 = b_comp
tt0 = tt
result = file_test(path0+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path0,my,mm,nanf,nz
endif else begin
 restore,path0+'dat/spectrum.sav'
endelse
restore,path1+'dat/'+strtrim(input)
b_comp1 = b_comp
tt1 = tt
endif
print,'restored!'
;#; modi da plottare
;if (d3 eq 1) then begin
;jin=69 ;usually corresponding to (1,-11)
;jfin=76 ;usually corresponding to (1,-5)
;endif else begin
;jin=1
;jfin=4
;endelse
;
filename_en=path0+'energy_compare_'+st0+'-'+st1+'.eps'
print,filename_en
phd_graphopener,filename=filename_en,xsize=20,ysize=15
meld = indgen(jfin-jin+1) * 0. + 1.
if (tmin ne 0.d) then begin
 bm0=  max(reform(b_comp0[tmin:idxmax0,jin:jfin]))
 bm1=  max(reform(b_comp1[tmin:idxmax1,jin:jfin]))
 bmm = max([bm0,bm1])
 yrange = [1.e-8,bmm]
endif else begin
 bm0 = max(reform(b_comp0[1:idxmax0,jin:jfin]))
 bm1 = max(reform(b_comp1[1:idxmax1,jin:jfin]))
 bmm = max([bm0,bm1])
 print,'AAAA',bm0,bm1
 yrange = [1.e-6,bmm]
endelse
print,yrange
xrange=[min(tt),tminmax]
print,'prova_xrange',xrange

m5=0.85
plot,tt0[tmin:idxmax0],b_comp0[tmin:idxmax0,jin] ,position=[m0,m2,m1,m5],/nodata,ytitle='W!DM!N',yr=yrange,ystyle=1,xr=xrange,xstyle=1,psym=0,/ylog

if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,tt0[tmin:idxmax0],b_comp0[tmin:idxmax0,q] ,col=colors0(q-jin),thick=th,psym=0
 oplot,tt1[tmin:idxmax1],b_comp1[tmin:idxmax1,q] ,col=colors0(q-jin),thick=th*2.,psym=0,linestyle=3
endfor
for q=nmds-1l,0l,-1l do begin
 col=colors(q)
; print,'q',q,q/2,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay
 xyouts,xbase+deltax*(q/2 + abs(q) mod 2),ybase-(q mod 2)*deltay,strings(q),col=col,/norm
endfor

endif else begin

for q=jin,jfin do begin
 print,'qqqq',q
 oplot,tt0[tmin:idxmax0],b_comp0[tmin:idxmax0,q],col=colors0(q),thick=th,psym=0
 oplot,tt1[tmin:idxmax1],b_comp1[tmin:idxmax1,q],col=colors0(q),thick=th*2.,psym=0,linestyle=3
endfor
xyouts,xbase,ybase+deltay,'SPC  r/a='+strtrim(string(float(raggior)/100.,format='(f6.3)')),col=0,/norm
xyouts,xbase+deltax,ybase+deltay,'h='+strtrim(string(-nanf(1),format='(i0)')),col=0,/norm
xyouts,xbase,ybase,'(1,h)',col=colors0(1),/norm
xyouts,xbase,ybase-deltay,'(2,2h)',col=colors0(2),/norm
xyouts,xbase+0.6*deltax,ybase,'(3,3h)',col=colors0(3),/norm
xyouts,xbase+0.6*deltax,ybase-deltay,'(4,4h)',col=colors0(4),/norm
endelse
xyouts,xbase+1.4*deltax,ybase-deltay,st0,col=0,/norm
xyouts,xbase+1.4*deltax,ybase-2*deltay,st1,col=colors1(1),/norm
xyouts,xbase+1.0*deltax,ybase-2*deltay,'dashed ',col=colors1(1),/norm

phd_graphcloser,filename_en

infiles = [filename,filename_v,filename_en]
infiles2 = strjoin(infiles, ' ')
outfile = path0+'compare_'+st0+'-'+st1+'.pdf'
spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
for i=0,n_elements(infiles) -1 do begin
 spawn, "rm -f "+infiles(i)
 print,infiles(i)
endfor
end
