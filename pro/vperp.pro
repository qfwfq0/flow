pro vperp,path,d3,itp0,itp1,deltaitp,nth,nzz,mymin,mymax

;#; programma che calcola le componenti perpendicolari e parallela della velocità, e anche il campo elettrico
;#; calcolo nello spazio reale
;path='/ricercatori/ft/specyl/veranda/flow/archive/rfp/3d/theta16/vario_dissipazione/varioH_n7_mp4.5e-3/'
print,path

itp0 = fix(itp0,type=3)
itp1 = fix(itp1,type=3)
nitp = itp1 - itp0 + 1
nsaved = ceil(float(nitp)/deltaitp)
aitp= indgen(nsaved)*deltaitp + itp0 
print,'itp',aitp

;#; Marco, mesh radiale della simulazione
mesh_size = read_mesh(path)
dx = 1.d/double(mesh_size)
pror1 = dindgen(mesh_size+1)/(mesh_size)
pror2 = dindgen(mesh_size+2)/(mesh_size) - 1.d/(2.*mesh_size)

;#; Marco, calcolo il profilo di resistività
dissipation = read_dissipation(path)
eta1 = make_array(mesh_size+1,/double,value=0.)
eta2 = make_array(mesh_size+2,/double,value=0.)
eta1 = dissipation.eta * (1. + (dissipation.alet-1)*pror1^dissipation.beet)^dissipation.gaet
eta2 = dissipation.eta * (1. + (dissipation.alet-1)*pror2^dissipation.beet)^dissipation.gaet

for k=0, n_elements(aitp)-1 do begin

 vperp_file=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/v_perp_'+strtrim(string(nth,'(i0)'))+'x'+strtrim(string(nzz,'(i0)'))+'_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'
 ;#; check whether v_perp file exists already
 existsav = file_info(vperp_file)
 if (existsav.exists eq 1) then begin
  continue
 endif
;#; ho bisogno della antiffft di B, J, v, da cui calcolo il campo elettrico E. Poi vperp e vpar sono semplici formulette dalla legge di Ohm.
 iitp = aitp(k)
 anti_fft_v=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/antifft_v_'+strtrim(string(n_elements(pror2),'(i0)'))+'x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'
 print,anti_fft_v
 anti_fft_dir=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/'
 call_procedure,'check_dir',anti_fft_dir
 anti_fft_b=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/antifft_b_'+strtrim(string(n_elements(pror2),'(i0)'))+'x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'
 anti_fft_j=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/antifft_j_'+strtrim(string(n_elements(pror2),'(i0)'))+'x'+strtrim(string(nth,'(i0)'))+'myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'

 ;#; check whether anti_fft_file exists
 hel = 1.
 deltaitp2 = 1
; mymin = 0
; mymax = 2
 exist = file_info(anti_fft_v)
 if (exist.exists eq 1) then begin
  restore,anti_fft_v
 endif else begin
  str = 'v'
  call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, iitp, iitp, deltaitp2, mymin, mymax
  restore,anti_fft_v
 endelse
 exist = file_info(anti_fft_b)
 if (exist.exists eq 1) then begin
  restore,anti_fft_b
 endif else begin
  str = 'b'
  call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, iitp, iitp, deltaitp2, mymin, mymax
  restore,anti_fft_b
 endelse
 exist = file_info(anti_fft_j)
 if (exist.exists eq 1) then begin
  restore,anti_fft_j
 endif else begin
  str = 'j'
  call_procedure, 'antifft_reconst_rfp', path, str, d3, nth, nzz, iitp, iitp, deltaitp2, mymin, mymax
  restore,anti_fft_j
 endelse
 
 help,vr
 help,jr
 help,jt
 help,jz
 ;#; eta in the 3d space
 eta3 = make_array(n_elements(vz[*,0,0]),n_elements(vz[0,*,0]),n_elements(vz[0,0,*]),/double,value=0.)
 
 er = vt*0.d
 et = vt*0.d
 ez = vt*0.d
 v_x_b = make_array(n_elements(vz[*,0,0]),n_elements(vz[0,*,0]),n_elements(vz[0,0,*]),3,/double,value=0.)
 j_x_b = make_array(n_elements(vz[*,0,0]),n_elements(vz[0,*,0]),n_elements(vz[0,0,*]),3,/double,value=0.)
 e_x_b = make_array(n_elements(vz[*,0,0]),n_elements(vz[0,*,0]),n_elements(vz[0,0,*]),3,/double,value=0.)
 
;#; Marco, interpolo certi vettori nella mesh x2
 vr2 = interpolate_x2_vec(reform(vr),path,d3) 
 br2 = interpolate_x2_vec(reform(br),path,d3) 
 jt2 = interpolate_x2_vec(reform(jt),path,d3) 
 jz2 = interpolate_x2_vec(reform(jz),path,d3) 

 ;#; modulus of B
 b2 = sqrt(br2^2.+bt^2.+bz^2.)
 help,jt2
 for p=0,nth-1 do begin
  if (p mod 50 eq 0) then begin
   print,'p=',p,'/',nth-1
  endif
  for q=0,nzz-1 do begin
    eta3[*,p,q] = eta2
    vxb = cross_product(vr2[*,p,q],vt[*,p,q],vz[*,p,q],br2[*,p,q],bt[*,p,q],bz[*,p,q])
    v_x_b[*,p,q,0] = vxb[*,0]
    v_x_b[*,p,q,1] = vxb[*,1]
    v_x_b[*,p,q,2] = vxb[*,2]
    jxb = cross_product(jr[*,p,q],jt2[*,p,q],jz2[*,p,q],br2[*,p,q],bt[*,p,q],bz[*,p,q])
    j_x_b[*,p,q,0] = jxb[*,0]
    j_x_b[*,p,q,1] = jxb[*,1]
    j_x_b[*,p,q,2] = jxb[*,2]
    er[*,p,q] = eta2 * jr[*,p,q] - v_x_b[*,p,q,0]
    et[*,p,q] = eta2 * jt2[*,p,q] - v_x_b[*,p,q,1]
    ez[*,p,q] = eta2 * jz2[*,p,q] - v_x_b[*,p,q,2]
    exb = cross_product(er[*,p,q],et[*,p,q],ez[*,p,q],br2[*,p,q],bt[*,p,q],bz[*,p,q])
    e_x_b[*,p,q,0] = exb[*,0]
    e_x_b[*,p,q,1] = exb[*,1]
    e_x_b[*,p,q,2] = exb[*,2]
  endfor
 endfor
 
 v_perp = make_array(n_elements(vz[*,0,0]),n_elements(vz[0,*,0]),n_elements(vz[0,0,*]),3,/double,value=0.)
 v_perp[*,*,*,0] = e_x_b[*,*,*,0]/b2 - eta3 * j_x_b[*,*,*,0]/b2
 v_perp[*,*,*,1] = e_x_b[*,*,*,1]/b2 - eta3 * j_x_b[*,*,*,1]/b2
 v_perp[*,*,*,2] = e_x_b[*,*,*,2]/b2 - eta3 * j_x_b[*,*,*,2]/b2
 
 save,filename=vperp_file,pror1,pror2,vecth,vecz,iitp,v_perp,er,et,ez
 
endfor ;#; end of cicle on ITP

;#; plot part
for k=0, n_elements(aitp)-1 do begin

  vperp_file=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/v_perp_'+strtrim(string(nth,'(i0)'))+'x'+strtrim(string(nzz,'(i0)'))+'_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.sav'
  restore,vperp_file

 v_perp2 = sqrt(v_perp[*,*,*,0]^2. + v_perp[*,*,*,1]^2. + v_perp[*,*,*,2]^2.)
 zed = 0
; mymin = 0
; mymax = 2

 letter="161B;"
 greektheta='!9' + String(letter) + '!X'
 letter2="160B;"
 greekpi='!9' + String(letter2) + '!X'
 letter="143B;"
 greekchi='!9' + String(letter) + '!X'
 letter="136B;"
 perpen0='!9' + String(letter) + '!X'
 letter="144B;"
 perpen1='!9' + String(letter) + '!X'
 yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
 ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]


 nlev = 100
 colors = findgen(nlev)/(nlev-1)*252.
 vperp2_pic=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/v_perp_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.eps'
 phd_graphopener,filename=vperp2_pic,xsize=12.,ysize=12.
  xrange=[-1.,1.]
  yrange = xrange
  contour,reform(v_perp2[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='v!D'+perpen0+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
  contour,reform(v_perp2[*,*,zed]),pror2#cos(vecth),pror2#sin(vecth),nlev=200,/fill,/overplot
  minv = min(v_perp2[*,*,zed])
  maxv = max(v_perp2[*,*,zed])
  xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
  xyouts,0.1,0.97,/norm, 'filled contour: v!D'+perpen0+'!N',chars=0.9
  xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(aitp(k),'(i0)')),chars=0.9
  xyouts,0.1,0.84,/norm, 'max v!D'+perpen0+'!N = '+strtrim(string(max(v_perp2[*,*,zed]),'(f6.4)')),chars=0.9
  colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)'
 phd_graphcloser,vperp2_pic

 vperpr_pic=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/v_perpr_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.eps'
 phd_graphopener,filename=vperpr_pic,xsize=12.,ysize=12.
  xrange=[-1.,1.]
  yrange = xrange
  contour,reform(v_perp[*,*,zed,0]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='v!D'+perpen0+',r!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
  contour,reform(v_perp[*,*,zed,0]),pror2#cos(vecth),pror2#sin(vecth),nlev=200,/fill,/overplot
  minv = min(v_perp2[*,*,zed,0])
  maxv = max(v_perp2[*,*,zed,0])
  xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
  xyouts,0.1,0.97,/norm, 'filled contour: v!D'+perpen0+',r!N',chars=0.9
  xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(aitp(k),'(i0)')),chars=0.9
  xyouts,0.1,0.84,/norm, 'max v!D'+perpen0+',r!N = '+strtrim(string(max(v_perp[*,*,zed,0]),'(f6.4)')),chars=0.9
  colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)'
 phd_graphcloser,vperpr_pic

 vperpt_pic=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/v_perpt_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.eps'
 phd_graphopener,filename=vperpt_pic,xsize=12.,ysize=12.
  xrange=[-1.,1.]
  yrange = xrange
  contour,reform(v_perp[*,*,zed,1]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='v!D'+perpen0+greektheta+'!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
  contour,reform(v_perp[*,*,zed,1]),pror2#cos(vecth),pror2#sin(vecth),nlev=200,/fill,/overplot
  minv = min(v_perp[*,*,zed,1])
  maxv = max(v_perp[*,*,zed,1])
  xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
  xyouts,0.1,0.97,/norm, 'filled contour: v!D'+perpen0+greektheta+'!N',chars=0.9
  xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(aitp(k),'(i0)')),chars=0.9
  xyouts,0.1,0.84,/norm, 'max v!D'+perpen0+greektheta+'r!N = '+strtrim(string(max(v_perp[*,*,zed,1]),'(f6.4)')),chars=0.9
  colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)'
 phd_graphcloser,vperpt_pic

 vperpz_pic=path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/v_perpz_myminmax'+strtrim(string(mymin,'(i0)'))+strtrim(string(mymax,'(i0)'))+'.eps'
 phd_graphopener,filename=vperpz_pic,xsize=12.,ysize=12.
  xrange=[-1.,1.]
  yrange = xrange
  contour,reform(v_perp[*,*,zed,2]),pror2#cos(vecth),pror2#sin(vecth),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,title='v!D'+perpen0+',z!N',xtitle='x/a',ytitle='y/a',/nodata,/iso,position=[0.19,0.15,0.82,0.82]
  contour,reform(v_perp[*,*,zed,2]),pror2#cos(vecth),pror2#sin(vecth),nlev=200,/fill,/overplot
  minv = min(v_perp[*,*,zed,2])
  maxv = max(v_perp[*,*,zed,2])
  xyouts,0.6,0.97,/norm, 'my!Dmin!N='+strtrim(string(mymin,format='(i0)'))+' my!Dmax!N='+strtrim(string(mymax,format='(i0)')),chars=0.9
  xyouts,0.1,0.97,/norm, 'filled contour: v!D'+perpen0+',z!N',chars=0.9
  xyouts,0.1,0.89,/norm, 'ITP:'+strtrim(string(aitp(k),'(i0)')),chars=0.9
  xyouts,0.1,0.84,/norm, 'max v!D'+perpen0+',z!N = '+strtrim(string(max(v_perp[*,*,zed,2]),'(f6.4)')),chars=0.9
  colorbar_marco,/vertical,color=250,c_colors=colors,min=minv,max=maxv,position=[0.83, 0.15, 0.86, 0.78],title='',/right,charsize=1.,divisions=10,format='(f7.4)'
 phd_graphcloser,vperpz_pic

 infiles = [vperp2_pic,vperpr_pic,vperpt_pic,vperpz_pic]
 infiles2 = strjoin(infiles, ' ')
 outfile = path+'dat/itp/'+strtrim(string(aitp(k),'(i0)'))+'/vperp_'+strtrim(string(mymin,format='(i0)'))+strtrim(string(mymax,format='(i0)'))+'.pdf'
 spawn, "convert -density 300 -antialias -rotate 0 -gravity center -background white "+infiles2+" -quality 100 " + outfile
 for i=0,n_elements(infiles) -1 do begin
; spawn, "rm -f "+infiles(i)
; print,infiles(i)
 endfor

endfor ;#; end of cicle on ITP

end   
