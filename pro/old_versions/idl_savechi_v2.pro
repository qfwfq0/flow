pro idl_savechi_v2,path,it,d3,tok,ieli,nth,nzz
;#; this routine computes the helical flux function.
;#; it can be computed on a 3D grids or on a 2D grid (both at th=cost,phi=cost)
;#; it computes Az=-integ(dr*Btheta) and r*Atheta=integ(dr*r*Bz)

 RR = read_major_radius(path)

;  restore,path+'dat/spectrum.sav'
 result = file_test(path+'dat/spectrum.sav')
 if (result ne 1) then begin
  call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
  if (m2d eq 0) then begin
   m2d = 1./ieli
  endif
 endif else begin
  restore,path+'dat/spectrum.sav'
 endelse
 print,'ieli=',ieli,' , it=',it,'  m2d=',m2d
;#; inputfile for vector potential
 input_a_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/acyl.sav'
 exx = file_test(input_a_file) 
 if (exx ne 1) then begin
  call_procedure, 'idl_savea',path,it,d3
  restore, input_a_file
 endif else begin
  restore, input_a_file
 endelse
 nr = n_elements(daz[*,0])
 pror = dindgen(nr) / (nr - 1)
 th = dindgen(nth) / (nth - 1) * 2. * !Dpi
 zz = dindgen(nzz) / (nzz - 1) * 2. * !Dpi * RR[0]
 rchi = make_array(nr,nth,nzz,/dcomplex,value=0.)
 rchi2 = make_array(nr,nth,nzz,/double,value=0.)
 rchi_ex = make_array(nr,nth,nzz,/dcomplex,value=0.)
 rchi_re = make_array(nr,nth,nzz,/double,value=0.)
 maxm = 5
 ii = dcomplex(0.d,1.d)
 dieli = double(ieli)

;#; begin 3d part
 if (tok eq 1) then begin
; if (d3 eq 1) then begin
  for z = 0 , nzz - 1 do begin
   if (z mod 32 eq 0) then print,'z',z
   for t = 0 , nth - 1 do begin
    rchi(*,t,z) =  daz[*,0] + dieli / RR[0] * dat[*,0]
    rchi2(*,t,z) =  rpdaz[*,0] + dieli / RR[0] * rpdat[*,0]
      for im = 1, maxm do begin
;#; I want just the nn=mm*ieli Fourier components of the vector potential.
       imm = im*m2d
       if (z eq 1 and t eq 1) then begin
        print,'m=',imm,' position in array=',im,' ',dieli
       endif
;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro
       fact =  2.*exp(ii * imm * ( th(t) + dieli *  zz(z) / RR[0]) )
       print,fact
       rchi(*,t,z) = rchi(*,t,z) + fact * ( daz[*,im] + dieli / RR[0] * dat[*,im] )

;#; calcolo usuale
;#; metto il fattore 1 invece che 2 perché ho calcolato il potenziale vettore rpda* e ipda* già con il fattore 2
       fcc =  1. * cos( imm * ( th(t) + dieli *  zz(z) / RR[0]) )
       fss =  - 1. * sin( imm * ( th(t) + dieli *  zz(z) / RR[0]) )
       rchi2(*,t,z) = rchi2(*,t,z) + fcc * (rpdaz[*,im] + dieli / RR[0] * rpdat[*,im] ) + fss * (ipdaz[*,im] + dieli / RR[0] * ipdat[*,im])
     endfor
   endfor
  endfor
; endif ;#;
  save_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel_v2.sav'
  save, filename=save_file,pror,th,zz,rchi,rchi2

 endif else begin;#; end tokamak part and begin RFP part

  maxm=2
  for z = 0 , nzz - 1 do begin
   if (z mod 32 eq 0) then print,'z',z,dieli
   for t = 0 , nth - 1 do begin
    rchi_ex(*,t,z) =  daz[*,0] - dieli / RR[0] * dat[*,0]
    rchi_re(*,t,z) =  rpdaz[*,0] - dieli / RR[0] * rpdat[*,0]
      for im = 1, maxm do begin
       inn =  im * dieli
       jmn = janz2(im,inn,my,mm,nanf,nzcon)
;       print,'check_RFP',im,inn,jmn,my,mm,nanf,nzcon
;#; I want just the nn=mm*ieli Fourier components of the vector potential.
       if (z eq 1 and t eq 1) then begin
        print,'m=',im,'n= ',inn,' position in array=',jmn,' ',dieli
       endif
;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro
       fact =  0.5*exp(ii * im * ( th(t) + dieli *  zz(z) / RR[0]) )
;       print,fact
       rchi_ex(*,t,z) = rchi_ex(*,t,z) + fact * ( daz[*,jmn] - dieli / RR[0] * dat[*,jmn] )

;#; calcolo usuale
       fcc =  0.5 * cos( im * ( th(t) + dieli *  zz(z) / RR[0]) )
       fss =  - 0.5 * sin( im * ( th(t) + dieli *  zz(z) / RR[0]) )
       rchi_re(*,t,z) = rchi_re(*,t,z) + fcc * (rpdaz[*,jmn] - dieli / RR[0] * rpdat[*,jmn] ) + fss * (ipdaz[*,jmn] - dieli / RR[0] * ipdat[*,jmn])
     endfor
   endfor
  endfor
  save_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chihel_v2.sav'
  save, filename=save_file,pror,th,zz,rchi_ex,rchi_re
 endelse

; save_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chi'+strtrim(string(dieli,format='(f3.2)'))+'.sav'
end
;
;   chit = rchi[*,*,0]
;   chiz = reform(rchi[*,0,*])
;
;   filename=path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/plot_chi.eps'
;   phd_graphopener,filename=filename,xsize=17.,ysize=5.
;
;   letter="161B;"
;   greektheta='!9' + String(letter) + '!X'
;   xr = [0.,1.]
;   yr = [0.,2.*!pi]
;   zr = yr*RR
;
;   !p.charsize=0.7
;   m1=0.06 & m2 =0.4 & m3=0.49 & m4=0.95
;   m5=0.17 & m6=0.90
;   xyouts,0.37,0.93,'helical flux function',/norm
;   !p.multi=[1,1,2,0,0]
;   contour,chit,pror,th,xrange=xr,yrange=yr,xtitle='r',ytitle=greektheta,title='z=0',position=[m1,m5,m2,m6],ystyle=5,yminor=1,/nodata
;   axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
;   axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
;   axis,yaxis=0,xrange=yr,col=0,  ytitle =greektheta,ystyle=1
;   axis,yaxis=1,xrange=yr,col=0,  ytitle ='',ychars=0.01,ystyle=1
;   contour,chit,pror,th,/overplot,nlevels=15
;
;   !p.multi=[1,1,2,0,0]
;   contour,transpose(chiz),zz,pror,/nodata,xrange=zr,yrange=xr,xtitle='z',ytitle='r',title=greektheta+'=0',position=[m3,m5,m4,m6],xstyle=5,xminor=1
;   axis,xaxis=0,xrange=zr,col=0,  xtitle ='z',xstyle=1
;   axis,xaxis=1,xrange=zr,col=0,  xtitle ='',xchars=0.01,xstyle=1
;   axis,yaxis=0,xrange=yr,col=0,  ytitle ='r',ystyle=1
;   axis,yaxis=1,xrange=yr,col=0,  ytitle ='',ychars=0.01,ystyle=1
;   contour,transpose(chiz),zz,pror,/overplot,nlev=12
;
;   phd_graphcloser,filename

