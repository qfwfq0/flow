pro idl_savechi,path,it,ieli,d3,nth,nzz
;#; this routine computes the helical flux function.
;#; it can be computed on a 3D grids or on a 2D grid (both at th=cost,phi=cost)
;#; it computes Az=-integ(dr*Btheta) and r*Atheta=integ(dr*r*Bz)

  ieli=fix(ieli)
  print,'ieli=',ieli,' , it=',it
;  if (d3 ne 1) then begin
;   d2 = 1
;  endif
;  path='./dat/'
;  it=260
  RR = read_major_radius(path)

;  restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
;  restore,path+'dat/itp/'+it+'/acyl.sav'
  restore,path+'/dat/ibprofiles.sav'

  nr = n_elements(pdbr[0,*,0])
  pror = findgen(nr) / (nr - 1)
  th = findgen(nth) / (nth - 1) * 2. * !Dpi
  zz = findgen(nzz) / (nzz - 1) * 2. * !Dpi * RR[0]
 
  rchi = make_array(nr,nth,nzz,/double,value=0.)
  maxm = 5

;#; compute the quantities needed for the equilibrium part of the helical flux function
  rpdaz_0 = make_array(n_elements(pdbr[0,*,0]),/double,value=0.)
  rpdat_0 = make_array(n_elements(pdbr[0,*,0]),/double,value=0.)
  rpdaz_0 = integral3(-reform(real_part(pdbt[it,*,0])),pror2)
  rpdat_0 = integral3(reform(real_part(pdbz[it,*,0])) * pror2,pror2)

;#; begin 3d part
  if (d3 eq 1) then begin
   for z = 0 , nzz - 1 do begin
    if (z mod 32 eq 0) then print,'z',z
    for t = 0 , nth - 1 do begin
;#; IMPORTANT: I put a minus sign because I'm integrating from edge to core
     rchi(*,t,z) = -rpdaz_0 + dieli / RR[0] * rpdat_0
      for im = 1, my do begin
;#; I want just the nn=mm*ieli Fourier components of the vector potential.
        nn = + im * ieli
        jmn = janz(im,nn,my,mm,nanf,nz)
        if (z eq 1 and t eq 1) then begin
         print,im,ieli*im,jmn
        endif
;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro!!
        fcos =  cos(im * th(t) + ieli * im * zz(z) / RR[0])
        fsin = - sin(im * th(t) + ieli * im * zz(z) / RR[0])
;        rchi(*,t,z) = rchi(*,t,z) + dcomplex(fcos * (im * rpdaz(*,jmn) - ieli / RR[0] * rpdat(*,jmn)) , fsin * (im * ipdaz(*,jmn) - ieli / RR[0] * ipdat(*,jmn)))
;        rchi(*,t,z) = rchi(*,t,z) + fcos * (im * rpdaz(*,jmn) - ieli / RR[0] * rpdat(*,jmn)) + fsin * (im * ipdaz(*,jmn) - ieli / RR[0] * ipdat(*,jmn))
;#; per i conti guardare il quaderno blu MHD14-15 chihel
        rchi(*,t,z) = rchi(*,t,z) + pror1 / dim * ( imaginary(reform(pdbr[it,*,jmn])) * fcos - (-real_part(reform(pdbr[it,*,jmn]))) * fsin)
      endfor
    endfor
   endfor
;#; end 3d part
  endif else begin

;#; begin 2d part
   dieli = double(ieli)

   for z = 0 , nzz - 1 do begin
    if (z mod 32 eq 0) then print,'z',zz(z)
    for t = 0 , nth - 1 do begin
;#; IMPORTANT: I put a minus sign because I'm integrating from edge to core
     rchi(*,t,z) = -rpdaz_0 + dieli / RR[0] * rpdat_0
;      for im = 1, my do begin
      for im = 1, maxm do begin
      dim = double(im)
;#; I want just the nn=mm*ieli Fourier components of the vector potential.
        nn = + im * ieli
;        print,'im',im,nn,nz
        jmn = janz2(im,nn,my,mm,nanf,nz)
        if (z eq 0 and t eq 0) then begin
         print,im,ieli*im,jmn
        endif
;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro!!
        fcos =  cos(dim * th(t) + dieli * dim * zz(z) / RR[0])
        fsin = - sin(dim * th(t) + dieli * dim * zz(z) / RR[0])
;#; per i conti guardare il quaderno blu MHD14-15 chihel
        rchi(*,t,z) = rchi(*,t,z) + pror1 / dim * ( imaginary(reform(pdbr[it,*,jmn])) * fcos - (-real_part(reform(pdbr[it,*,jmn]))) * fsin)

;#; old methods with problems.
;        rchi(*,t,z) = rchi(*,t,z) + dcomplex(fcos * (im * rpdaz(*,jmn) - ieli / RR[0] * rpdat(*,jmn)) , fsin * (im * ipdaz(*,jmn) - ieli / RR[0] * ipdat(*,jmn)))
;       rchi(*,t,z) = rchi(*,t,z) + fcos * (dim * rpdaz(*,jmn) - dieli / RR[0] * rpdat(*,jmn)) + fsin * (dim * ipdaz(*,jmn) - dieli / RR[0] * ipdat(*,jmn))
      endfor
    endfor
   endfor


  endelse
;#; end 2d part

   save_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/chi'+strtrim(string(ieli,format='(i0)'))+'.sav'
   save, filename=save_file,pror,th,zz,rchi

   chit = rchi[*,*,0]
   chiz = reform(rchi[*,0,*])

   filename=path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/plot_chi.eps'
   phd_graphopener,filename=filename,xsize=17.,ysize=5.

   letter="161B;"
   greektheta='!9' + String(letter) + '!X'
   xr = [0.,1.]
   yr = [0.,2.*!pi]
   zr = yr*RR

   !p.charsize=0.7
   m1=0.06 & m2 =0.4 & m3=0.49 & m4=0.95
   m5=0.17 & m6=0.90
   xyouts,0.37,0.93,'helical flux function',/norm
   !p.multi=[1,1,2,0,0]
   contour,chit,pror,th,xrange=xr,yrange=yr,xtitle='r',ytitle=greektheta,title='z=0',position=[m1,m5,m2,m6],ystyle=5,yminor=1,/nodata
   axis,xaxis=0,xrange=xr,col=0,  xtitle ='r'
   axis,xaxis=1,xrange=xr,col=0,  xtitle ='',xchars=0.01
   axis,yaxis=0,xrange=yr,col=0,  ytitle =greektheta,ystyle=1
   axis,yaxis=1,xrange=yr,col=0,  ytitle ='',ychars=0.01,ystyle=1
   contour,chit,pror,th,/overplot,nlevels=15

   !p.multi=[1,1,2,0,0]
   contour,transpose(chiz),zz,pror,/nodata,xrange=zr,yrange=xr,xtitle='z',ytitle='r',title=greektheta+'=0',position=[m3,m5,m4,m6],xstyle=5,xminor=1
   axis,xaxis=0,xrange=zr,col=0,  xtitle ='z',xstyle=1
   axis,xaxis=1,xrange=zr,col=0,  xtitle ='',xchars=0.01,xstyle=1
   axis,yaxis=0,xrange=yr,col=0,  ytitle ='r',ystyle=1
   axis,yaxis=1,xrange=yr,col=0,  ytitle ='',ychars=0.01,ystyle=1
   contour,transpose(chiz),zz,pror,/overplot,nlev=12

   phd_graphcloser,filename

end
