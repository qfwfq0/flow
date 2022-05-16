pro show_time_chihel , path,itp0,itp1,inteli,deltaitp,itp

loadct,'0',file='~/idl/colors.tbl'
outfile = path+'dat/itp/chitot_ieli'+strtrim(string(inteli,format='(i0)'))+'_'+strtrim(string(itp0,format='(i0)'))+'-'+strtrim(string(itp1,format='(i0)'))+'-'+strtrim(string(deltaitp,format='(i0)'))+'.sav'
restore,outfile

phi = 0
window,1
nlevs=45
for q = 0, n_elements(rchi_tot[*,0,0,0]) - 1 do begin
 contour,reform(rchi_tot[q,*,*,phi]),pror1#cos(th),pror1#sin(th),/iso,nlev=nlevs,color=255,/fill,xtitle='x',ytitle='y',title='chi hel, time='+strtrim(string(itp(q),'(f9.4)'))
 contour,reform(rchi_tot[q,*,*,phi]),pror1#cos(th),pror1#sin(th),nlev=nlevs,c_color=255,/overplot
 wait,1.
endfor

end
