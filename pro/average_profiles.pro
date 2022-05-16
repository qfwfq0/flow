pro average_profiles, path,d3

input='imf_bprofiles.sav'
restore,path+'dat/'+strtrim(input)
;restore,path+'dat/spectrum.sav'
print,path
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse

;#; Marco, plot bz(a) to detect [tmin,tmax] over which average the profiles.
;#; modi da plottare
if (d3 eq 1) then begin
jin=69 ;usually corresponding to (1,-11)
jfin=76 ;usually corresponding to (1,-5)
endif else begin
jin=1
jfin=2
endelse

rds = 0.88
raggior = min(where(pror1 gt rds))
raggiot = min(where(pror2 gt rds))
colors=[180,0,250,20,115,75,220,50]
th=3.

;#; scelgo di fare le medie solo su tempi equispaziati!
 deltat = tt-shift(tt,1)
 deltat=deltat(where(deltat gt 0.))
 delta = max(deltat)
 pick = where( ( abs(tt mod delta[0]) lt 1.e-5 ) or ( abs(tt mod delta[0] - delta[0]) lt 1.e-5) )

;#; ridefinisco i vettori dei campi
 br = br[pick,*,*,*]
 bt = bt[pick,*,*,*]
 bz = bz[pick,*,*,*]
 tt = tt[pick]

set_plot,'x'
window,0
plot,bz[*,raggiot,73,0],xtitle='time steps',ytitle='bz(a)',thick=0.01,/nodata
if (d3 eq 1) then begin
for q=jin,jfin do begin
 oplot,bz[*,raggiot,q,0] * 2.,col=colors(q-jin),thick=th
endfor
endif else begin
for q=jin,jfin do begin
 oplot,bz[*,raggiot,q,0] * 2.,col=colors(q),thick=th
endfor
endelse

tmin=''
tmax=''
read,tmin,prompt='write timestep minimum:  '
read,tmax,prompt='write timestep maximum:  '
print,'tmin=',tmin,'  tmax=',tmax

wdelete,0

;#; Marco, creo array con i tre profili mediati di campo magnetico

av_br = make_array(n_elements(br[0,*,0,0]),n_elements(br[0,0,*,0]),/double,value=0.d)
av_br2 = make_array(n_elements(br[0,*,0,0]),n_elements(br[0,0,*,0]),/double,value=0.d)
av_bt = make_array(n_elements(bt[0,*,0,0]),n_elements(bt[0,0,*,0]),/double,value=0.d)
av_bz = make_array(n_elements(bz[0,*,0,0]),n_elements(bz[0,0,*,0]),/double,value=0.d)

tmin=fix(tmin)
tmax=fix(tmax)
if (tmax ge n_elements(br[*,0,0,0])) then begin
 tmax = n_elements(br[*,0,0,0]) - 1
endif
if (tmin ge n_elements(br[*,0,0,0])) then begin
 read,tmin,prompt='TRY AGAIN! write timestep minimum:  '
 tmin=fix(tmin)
endif

delta = tmax-tmin
print,delta
help,br[tmin:tmax,0,0,0]
for p=0, n_elements(br[0,0,*,0])-1 do begin
 for q=0, n_elements(br[0,*,0,0])-1 do begin
  av_br[q,p] = total(abs(br[tmin:tmax,q,p,0]) /bt[tmin:tmax,101,0,0] ) / delta
 endfor
 for q=0, n_elements(bt[0,*,0,0])-1 do begin
  av_bt[q,p] = total(abs(bt[tmin:tmax,q,p,0]) /bt[tmin:tmax,101,0,0] ) / delta
  av_bz[q,p] = total(abs(bz[tmin:tmax,q,p,0]) /bt[tmin:tmax,101,0,0] ) / delta
 endfor
endfor

savfile=path+'dat/av_bfield.sav'
save,filename=savfile,pror1,pror2,tmin,tmax,av_br,av_bt,av_bz

end
