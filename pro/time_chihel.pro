function integral2,vect,xx
;#; general purpose function to calculate the integral of a function (vect) along the mesh x2 (vect)
;#; x2=findgen(lx+2)/lx-1./(2.*lx)
;#; the result will be stored in the mesh x1(i) = findgen(lx+1)/lx

 integ = make_array(n_elements(xx)-1,n_elements(vect(0,*)),/double,value=0.)
 yy = findgen(n_elements(xx)-1) / (n_elements(xx)-2)
 dx = 1.d / (n_elements(xx)-2)
 
 for q=0, n_elements(vect(0,*))-1 do begin
  integ(0,q) = 0.d
  for i=1,n_elements(yy)-1 do begin
   integ(i,q) = integ(i-1,q) + vect(i,q) * dx
  endfor
 endfor
 return,integ
end

pro time_chihel, path, d3, tok, hel, itp0, itp1, deltaitp

RR = read_major_radius(path)
if (RR(0) lt 1e-1) then begin
 RR(0)=4.d
endif
;path = './'
;d3 = 0
;hel = 10.
;#; understand which kind of energy is being computed'
; restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
 print,nz
 nzcon=nz
 nzcon[0] = nz[0] + 1
 i2d=1 
 m2d=2
 n2d=1
endelse
 mf_field='dat/ibprofiles.sav'
 dum = file_test(path+mf_field)
 print,'cacca',path+mf_field
 if (dum ne 1) then begin
  call_procedure, 'isavemodeprofiles',path,0,d3
 endif
 cmpt= test_up_to_date2(path,mf_field)
 print,'cmpt',cmpt
 if (cmpt gt 0) then begin
  call_procedure, 'isavemodeprofiles',path,0,d3
 endif
 
!p.charsize=1.75
 restore, mf_field
; print, 'restored!'
; print, 'helicity=', hel
 hel = double(hel)

 nr1 = n_elements(pror1)
 nr2 = n_elements(pror2)
 nth = 256
 vectth = dindgen(nth)/(nth-1) * 2. * !Dpi
 ;vecz = dindgen(nz)/(nz-1) * 2. * !Dpi * RR[0],;,
 
 outfile = path+'/dat/time_chihel.sav'
 if (itp1 eq -1) then begin 
  itp1=n_elements(pdbr[*,0,0,0])-1
  itp0=n_elements(pdbr[*,0,0,0])-1
  deltaitp=1
 endif
 if (itp0 eq itp1) then begin
  outfile = path+'dat/itp/'+strtrim(string(itp0,format='(i0)'))+'/chihel.sav'
  print,'need to create the directory for helical flux function file for itp=',itp0
  dir =  path+'dat/itp/'+strtrim(string(itp0,format='(i0)'))+'/'
  print,'dir',dir
  call_procedure,'check_dir',dir 
  ex_file = file_test(outfile)
 endif


;#; if cycle on the file time_chihel.sav existence
 dum = file_test(outfile)
 tme_sav = file_info(outfile)
 time_sav = tme_sav.mtime 
 tme_profsav = file_info(path+mf_field)
 time_profsav = tme_profsav.mtime 
; if (dum ne 1 or time_sav lt time_profsav) then begin
 print,'2d reconstruction with hel= ',hel
 if (tok eq 1) then begin
  maxm = 7
 endif else begin
 print,'rfp part'
 if (d3 eq 1) then begin
  maxm = 2
 endif else begin
  maxm=5
 endelse
 endelse

;#; Marco, deal with time
 if (itp0 eq !null) then begin
  itp0 = 0
 endif
 if (itp1 eq !null) then begin
  itp1 = n_elements(pdt) - 1
 endif
 
 if (itp0 eq !null) then begin
  print,'define itp0!'
 endif
 itp0 = fix(itp0,type=3)
 itp1 = fix(itp1,type=3)
 nitp = itp1 - itp0 + 1
 nsaved = ceil(float(nitp)/deltaitp)
 itp= indgen(nsaved)*deltaitp + itp0 

; rchi = make_array(n_elements(pdbz[*,0,0]),nr1,nth,/double,value=0.)
 rchi2 = make_array(nsaved,nr1,nth,/double,value=0.)
 help,rchi2
   letter="161B;"
   greektheta='!9' + String(letter) + '!X'
   letter2="160B;"
   greekpi='!9' + String(letter2) + '!X'
   letter="143B;"
   greekchi='!9' + String(letter) + '!X'
; print,'helicity=',hel
;#; cycle on time to ceate vector potential

 if (tok eq 1) then begin
 if (i2d eq 1) then begin
;  m2d=2
  for q = 0, nsaved - 1 do begin
;#; compute the quantities needed for the equilibrium part of the helical flux function
  rpdaz_0 = make_array(nsaved,n_elements(pdbr[0,*,0]),/double,value=0.)
  rpdat_0 = make_array(nsaved,n_elements(pdbr[0,*,0]),/double,value=0.)
  rpdaz_0[q,*] = integral3(-reform(real_part(pdbt[itp(q),*,0])),pror2)
  rpdat_0[q,*] = integral3(reform(real_part(pdbz[itp(q),*,0])) * pror2,pror2)
  dum = make_array(nr1,nth,/double,value=0.)
  dum = dum * 0.d

    for t = 0 , nth - 1 do begin
;#; chissà perché per il tokamak ci va un segno diverso rispetto a rfp. indagare.
;#; Marco, 07/04/2021 versione originale
;     dum(*,t) = -rpdaz_0[q,*] - hel / RR[0] * rpdat_0[q,*]
;#; Marco, 07/04/2021 prove
     dum(*,t) = rpdaz_0[q,*] - hel / RR[0] * rpdat_0[q,*]
  ;   dum(*,t) = reform(dmrpdaz(*,0)) - hel / RR[0] * reform(dmrpdat(*,0))
      for im = 1, maxm do begin
;#; I want just the nn=mm*hel Fourier components of the vector potential.
;#; Marco, 2D case!
        mmm = mm(im)*m2d
        nn = double(- mmm * hel)
      
        if (im mod m2d eq 0) then begin
         nn = double(nanf(im))
        jmn3 = janz3(im,nn,my,mm,nanf,nz)
        jmn2 = janz2(im,nn,my,mm,nanf,nz)
        jmn = janz(im,nn,my,mm,nanf,nz)
        jmn3 = im/m2d
;        print,'tok2d spectrum',im,nn,jmn3,jmn2,jmn
        dim = double(im)
         if ( t eq 0) then begin
          print,'tok2d m=',im,nn,itp(q),jmn3,jmn2,jmn
         endif
;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro!!
        fcos =  cos(im * vectth)
        fsin = - sin(im * vectth)
;#; per i conti guardare il quaderno blu MHD14-15 chihel
        dum(*,t) = dum(*,t) + pror1 / dim * ( imaginary(reform(pdbr[itp(q),*,jmn3])) * fcos(t) - (-real_part(reform(pdbr[itp(q),*,jmn3]))) * fsin(t))
       
        endif else begin
        endelse
      endfor
    endfor
   rchi2(q,*,*) = dum
   print,'BBB',itp(q)
   rchi = reform(rchi2(q,*,*))
   dir =  path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/'
   call_procedure,'check_dir',dir 
   otfile = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel.sav'
   save,filename=otfile,rchi,pror1,vectth,pdt,itp
   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel.eps'
   letter="161B;"
   greektheta='!9' + String(letter) + '!X'
   letter2="160B;"
   greekpi='!9' + String(letter2) + '!X'
   letter="143B;"
   greekchi='!9' + String(letter) + '!X'
   yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
   ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

;#; find separatrix
   aa = max(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   cc = min(rchi[bb[0],*]) ;#; separatrix!
   print,'sep_lev=',cc
   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
   nisla=5
   dd=findgen(nisla)/(nisla-1)*abs((cc-aa))+cc
   dd = dd(uniq(dd))
   print,'dd',dd
;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;   print,'ee',ee
   phd_graphopener,filename=image_file,xsize=15.,ysize=10.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange=[0.,1.]
    yrange=[0.,2.*!pi]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1,vectth,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=200,col=0,yticks=5,ytickv=yt0,ytickname=ytpi,title=greekchi+'!Dhel!N',xtitle='r/a',ytitle=greektheta
    contour,reform(rchi[*,*]),pror1,vectth,lev=cc[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,levels=dd,col=253,thick=1.,/overplot

   phd_graphcloser,image_file
   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel_xy.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/iso,position=[0.27,0.21,0.94,0.94]
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=214,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
   contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=253,thick=1.,/overplot

   phd_graphcloser,image_file

   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel_xy_fill.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,rchi[*,*],pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/iso,/fill,position=[0.15,0.15,0.90,0.90]
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=253,thick=1.,/overplot,/fill
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot

   phd_graphcloser,image_file
   endfor
   outfile = path+'/dat/time_chihel_'+strtrim(string(itp0,'(i0)'))+strtrim(string(itp1,'(i0)'))+'.sav'
   save,filename=outfile,rchi2,pror1,vectth,pdt,itp

  endif else begin;#; end of 2D TOK case
   print,'code the TOK3D case'
   for q = 0, nsaved - 1 do begin
; #; compute the quantities needed for the equilibrium part of the helical flux function
   rpdaz_0 = make_array(nsaved,n_elements(pdbr[0,*,0]),/double,value=0.)
   rpdat_0 = make_array(nsaved,n_elements(pdbr[0,*,0]),/double,value=0.)
   rpdaz_0[q,*] = integral3(-reform(real_part(pdbt[itp(q),*,0])),pror2)
   rpdat_0[q,*] = integral3(reform(real_part(pdbz[itp(q),*,0])) * pror2,pror2)
   dum = make_array(nr1,nth,/double,value=0.)
   dum = dum * 0.d
    for t = 0 , nth - 1 do begin
; #; chissà perché per il tokamak ci va un segno diverso rispetto a rfp. indagare.
;#; Marco, 01/03/2021 versione originale
;     dum(*,t) = -rpdaz_0[q,*] - hel / RR[0] * rpdat_0[q,*]
;#; Marco, 01/03/2021 prove
     dum(*,t) = rpdaz_0[q,*] - hel / RR[0] * rpdat_0[q,*]
   ;   dum(*,t) = reform(dmrpdaz(*,0)) - hel / RR[0] * reform(dmrpdat(*,0))
     for im = 1, maxm do begin
; #; I want just the nn=mm*hel Fourier components of the vector potential.
; #; Marco, 2D case!
;         mmm = mm(im)*m2d
;         print,mmm
;         nn = double(- mmm * hel)
;       if ((im mod 2) eq 0) then begin 
        nn = im*hel
        jmn2 = janz2(im,nn,my,mm,nanf,nzcon)
;        print,'tok3d spectrum',hel,im,nn,jmn2
        dim = double(im)
       if ( t eq 0) then begin
        print,'tok3d m=',im,nn,itp(q),jmn2
       endif
; #; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro!!
        fcos =  cos(im * vectth)
        fsin = - sin(im * vectth)
; #; per i conti guardare il quaderno blu MHD14-15 chihel
        dum(*,t) = dum(*,t) + pror1 / dim * ( imaginary(reform(pdbr[itp(q),*,jmn2])) * fcos(t) + real_part(reform(pdbr[itp(q),*,jmn2])) * fsin(t) )
;       endif
     endfor
    endfor
    rchi2(q,*,*) = dum
    rchi = reform(rchi2(q,*,*))
    dirdd = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/'
    call_procedure,'check_dir',dirdd
    otfile = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel.sav'
    save,filename=otfile,rchi,pror1,vectth,pdt,itp

;#; find separatrix
   aa = max(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   cc = min(rchi[bb[0],*]) ;#; separatrix!
   print,'sep_lev=',cc
   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
   nisla=5
   dd=findgen(nisla)/(nisla-1)*abs((cc-aa))+cc
   dd = dd(uniq(dd))
   print,'dd',dd
   ytickv=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
   ytickname=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel_xy.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
    xrange2=[-1.,1.]
    contour,rchi[*,*],pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/iso,position=[0.15,0.15,0.90,0.90]
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=53,thick=1.,/overplot

   phd_graphcloser,image_file

   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel_xy_fill.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
    xrange2=[-1.,1.]
    contour,rchi[*,*],pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/iso,/fill,position=[0.15,0.15,0.90,0.90]
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=cc[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=53,thick=1.,/overplot,/fill
   phd_graphcloser,image_file

   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel_fill.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=7.
    xrange3=[0.,1.]
    yrange3=[0.,2.*!PI]
    contour,rchi[*,*],pror1,vectth,xrange=xrange3,yrange=yrange3,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N',xtitle='r/a',ytitle=greektheta,/fill,position=[0.15,0.15,0.90,0.88]
;    contour,reform(rchi[*,*]),pror1,vectth,lev=cc[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,levels=dd,col=53,thick=1.,/overplot,/fill
   phd_graphcloser,image_file

   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=7.
    xrange3=[0.,1.]
    yrange3=[0.,2.*!pi]
;    print,'abc',yrange3
;    print,ytickv
;    print,ytickname
    plot,[0.],xrange=xrange3,yrange=yrange3,xstyle=1,ystyle=1,title=greekchi+'!Dhel!N',xtitle='r/a',ytitle=greektheta,position=[0.15,0.15,0.90,0.88],/nodata
    contour,rchi[*,*],pror1,vectth,nlev=20,col=0,/overplot
;    contour,reform(rchi[*,*]),pror1,vectth,lev=cc[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,levels=dd,col=53,thick=1.,/overplot
   phd_graphcloser,image_file

   endfor ;#; end cycle on q

  endelse;#; end of 3D TOK case

endif else begin ;#; RFP case

  print,'RFP_case'
  for q = 0, nsaved - 1 do begin
  otfile = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel.sav'
;  print,'need to create the directory for helical flux function file for itp=',itp(q)
  dir =  path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/'
;  print,'dir',dir
  call_procedure,'check_dir',dir 
  ex_file = file_test(otfile)
;#; compute the quantities needed for the equilibrium part of the helical flux function
  rpdaz_0 = make_array(nsaved,n_elements(pdbr[0,*,0]),/double,value=0.)
  rpdat_0 = make_array(nsaved,n_elements(pdbr[0,*,0]),/double,value=0.)
  rpdaz_0[q,*] = integral3(-reform(real_part(pdbt[itp(q),*,0])),pror2)
  rpdat_0[q,*] = integral3(reform(real_part(pdbz[itp(q),*,0])) * pror2,pror2)
  dum = make_array(nr1,nth,/double,value=0.)
  dum = dum * 0.d
    for t = 0 , nth - 1 do begin
;#; IMPORTANT: I put a minus sign because I'm integrating from edge to core
     dum(*,t) = -rpdaz_0[q,*] + hel / RR[0] * rpdat_0[q,*]
      for im = 1, maxm do begin
;#; I want just the nn=mm*hel Fourier components of the vector potential.
;#; Marco, 2D case!
       if (d3 eq 0) then begin
        mmm = mm(im)
        jmn = im
       endif else begin
        mmm = mm(im)
        inn = mmm * hel
        jmn = janz2(mmm,inn,my,mm,nanf,nzcon)
;#; Marco, 3D case!
       endelse
        dim = double(mmm)
        if (q eq 0 and t eq 0) then begin
         print,'chihel modes',im,mmm,hel*mmm,jmn,hel,maxm
        endif
;#; Marco, tolgo il fattore due perché è già presente in isavemodeprofiles.pro!!
        fcos =  cos(mmm * vectth)
        fsin = - sin(mmm * vectth)
;#; per i conti guardare il quaderno blu MHD14-15 chihel
        dum(*,t) = dum(*,t) + pror1 / dim * ( imaginary(reform(pdbr[itp(q),*,jmn])) * fcos(t) - (-real_part(reform(pdbr[itp(q),*,jmn]))) * fsin(t))
       
      endfor
    endfor
    rchi2(q,*,*) = dum
    rchi = reform(rchi2(q,*,*))
; print,'outfile_rfp',outfile
   save,filename=otfile,rchi,pror1,vectth,pdt,itp
;   print,'image_file',image_file
   letter="161B;"
   greektheta='!9' + String(letter) + '!X'
   letter2="160B;"
   greekpi='!9' + String(letter2) + '!X'
   letter="143B;"
   greekchi='!9' + String(letter) + '!X'
   yt0=[0.,!pi/2.,!pi,3.*!pi/2.,2.*!pi]
   ytpi=['0',greekpi+'/2',greekpi,'3'+greekpi+'/2','2'+greekpi]

;#; find separatrix
;   q=0
   aa = max(rchi[*,*],idx) ;#; O-point
   bb = array_indices(rchi[*,*],idx)
   ee0 = min(rchi[bb[0],*],eeidx)
   ee = max(rchi[*,eeidx])
   print,bb
;#; altra tecnica per trovare la separatrice guardando il minimo della derivata e poi cercando il massimo della chi a sinistra del minimo della derivata
;   if (bb[0] eq 100) then begin
;    drchi = rchi*0.d
;    for k=0,n_elements(rchi[0,*])-1 do begin
;     drchi[*,k] = deriv(pror1,rchi[*,k])
;    endfor
;     aa=min(drchi,idx)
;     bb=array_indices(drchi,idx)
;   ee0 = min(rchi[bb[0],*],eeidx)
;   ee = max(rchi[*,eeidx])
;   kk=max(rchi[0:bb[0],bb[1]],kidx)
;   endif else begin
;   kk=1.
;   endelse
;
;   print,'sep_lev=',ee,kk
;   print,'Opoint_lev=',aa;,min(cc[0],aa[0])
;   nisla=5
;   dd=findgen(nisla)/(nisla-1)*abs((ee-aa))+ee
;   print,'dd',dd


   if (bb[0] eq 100) then begin
    drchi = rchi*0.d
    for k=0,n_elements(rchi[0,*])-1 do begin
     drchi[*,k] = deriv(pror1,rchi[*,k])
    endfor
     aa=min(drchi,idx)
     bb=array_indices(drchi,idx)
   ee0 = min(rchi[bb[0],*],eeidx)
   ropoint = min(rchi[bb[0]:*,bb[1]],idxopoint)
   ee = max(rchi[*,eeidx])
   kk=max(rchi[0:bb[0],bb[1]],kidx)
   endif else begin
   kk=1.
   idxopoint=0
   endelse


   nisla=10
   dd=findgen(nisla)/(nisla-1)*abs((ee-aa))+ee
   print,'dd',dd
   alfa=0.315
   dd1=findgen(nisla)/(nisla-1)*abs((ee-aa))*(1+alfa)+ee-alfa/2.*abs(ee-aa)
   dd1 = dd1(uniq(dd1))
   print,'dd1',dd1
   print,'ee',ee,'aa',aa
   print,pror1(bb(0)),vectth(eeidx)
   print,'opoint',pror1(bb[0]+idxopoint(0)),vectth(bb[1])

;   ee=findgen(10)/9*abs((cc-aa))+min(cc,aa)
;   print,'ee',ee
   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel.eps'
   phd_graphopener,filename=image_file,xsize=15.,ysize=10.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange=[0.,1.]
    yrange=[0.,2.*!pi]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1,vectth,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,col=0,yticks=5,ytickv=yt0,ytickname=ytpi,title=greekchi+'!Dhel!N',xtitle='r/a',ytitle=greektheta
    contour,reform(rchi[*,*]),pror1,vectth,nlev=20,col=0,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,lev=ee[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1,vectth,levels=dd,col=23,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1,vectth,levels=kk,col=43,thick=1.,/overplot
;    contour,shift(rchi[*,*],0,66),pror1,vectth,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,nlev=20,col=0,yticks=5,ytickv=yt0,ytickname=ytpi,title=greekchi+'!Dhel!N',xtitle='r/a',ytitle=greektheta
;    contour,shift(rchi[*,*],0,66),pror1,vectth,nlev=20,col=0,thick=1.,/overplot
;    contour,shift(rchi[*,*],0,66),pror1,vectth,lev=ee[0],col=254,thick=3.,/overplot
;    contour,shift(rchi[*,*],0,66),pror1,vectth,lev=aa[0]*0.99999,col=253,thick=3.,/overplot
;;    contour(shift(rchi[*,*],0,66)),pror1,vectth,levels=dd,col=23,thick=1.,/overplot
;    contour,shift(rchi[*,*],0,66),pror1,vectth,levels=kk,col=43,thick=1.,/overplot

   phd_graphcloser,image_file
   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel_xy.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
;#; specifico per un caso paper riconnessione 2020 / simulazione 2D h9 s100m1000 qinfty momsour 0
;   rchi=shift(rchi,0,-43) 
;   rchi=shift(rchi,0,-40) 
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=5,ystyle=5,nlev=20,col=0,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/iso
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=254,thick=3.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=aa[0]*0.99999,col=253,thick=3.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=223,thick=1.,/overplot
    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot

   phd_graphcloser,image_file

   image_file = path+'dat/itp/'+strtrim(string(itp(q),format='(i0)'))+'/chihel_xy_fill.eps'
   phd_graphopener,filename=image_file,xsize=12.,ysize=12.
;    contour,reform(rchi[q,*,*]),pror1#cos(vectth),pror1#sin(vectth),/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
;    contour,reform(rchi[q,*,*]),pror1,vectth,/iso,nlev=50,color=col(q),/fill,title=string(q),/nodata
    xrange2=[-1.,1.]
    xrange2=[-0.75,0.75]
;    plot,[0],xrange=xrange,yrange=yrange,/nodata,xtitle='r/a',ytitle=greektheta,ystyle=1,xstyle=1
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N itp='+strtrim(string(itp(q),format='(f5.2)')),xtitle='x/a',ytitle='y/a',/iso,/fill,position=[0.24,0.15,0.94,0.94],ycharsize=0.01,xticks=4,xtickv=[-0.5,0.,0.5,1.],xtickname=['-0.5',' ','0.5','1']
  help,rchi
;#; specifico per un caso paper riconnessione 2020 / simulazione 2D /ricercatori/ft/specyl/veranda/flow/archive_flow/rfp/2d/theta16/momsour/h9/s100m1000_Qinfty_more_harmonics/
;   rchi=shift(rchi,0,-63) 

nlev=40
c_color=findgen(nlev)/(nlev-1)*254.
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=20,col=0,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/iso,/fill,position=[0.21,0.18,0.95,0.95],ycharsize=0.01,xticks=4,xtickv=[-0.5,0.,0.5,1.],xtickname=['-0.5',' ','0.5','1'],chars=1.7
   alfa=1
   contour,reform(rchi[*,0:-alfa]),pror1#cos(vectth[0:-alfa]),pror1#sin(vectth[0:-alfa]),xrange=xrange2,yrange=xrange2,xstyle=1,ystyle=1,nlev=nlev,c_col=c_color,c_thick=3.5,title=greekchi+'!Dhel!N',xtitle='x/a',ytitle='y/a',/iso,position=[0.21,0.18,0.95,0.95],chars=1.7;,xticks=4,xtickv=[-0.5,0.,0.5,1.],xtickname=['-0.5',' ','0.5','1']
;loadct,7,/silent
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),lev=ee[0],col=25,thick=2.,/overplot
    contour,reform(rchi[*,0:-alfa]),pror1#cos(vectth[0:-alfa]),pror1#sin(vectth[0:-alfa]),lev=aa[0]*0.99999,col=25,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd[3],col=25,thick=1.,/overplot
    dd_col=findgen(n_elements(dd1))/(n_elements(dd)-1)*50+200.
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=dd,col=25,thick=1.,/overplot
    contour,reform(rchi[*,0:-alfa]),pror1#cos(vectth[0:-alfa]),pror1#sin(vectth[0:-alfa]),levels=dd1,c_col=dd_col,c_thick=2.5,thick=1.,/overplot
;    contour,reform(rchi[*,*]),pror1#cos(vectth),pror1#sin(vectth),levels=kk,col=43,thick=1.,/overplot
loadct,'0',file='~/idl/colors.tbl',/silent

   phd_graphcloser,image_file

   endfor
   outfile = path+'/dat/time_chihel_'+strtrim(string(itp0,'(i0)'))+strtrim(string(itp1,'(i0)'))+'.sav'
   save,filename=outfile,rchi2,pror1,vectth,pdt,itp
  endelse

;#; save
; endif else begin
;#; restore already computed
; restore,filename=outfile
; endelse


 
;  call_procedure, 'show_time_chihel',path,itp0,itp1
;pro show_time_chihel , path,itp0,itp1,inteli,deltaitp,itp
  

;loadct,'0',file='~/idl/colors.tbl'
;col=findgen(n_elements(rchi[*,0,0]))/(n_elements(rchi[*,0,0]) - 1)*253.
;;window,1
;; for q = 0, n_elements(br[*,0,0,0]) - 1 do begin
;;    contour,reform(bbr[q,*,*,0]),pror1,vectth,nlev=20,color=col(q)
;;   wait,0.02
;; endfor
;;window,2
;; for q = 0, n_elements(br[*,0,0,0]) - 1 do begin
;;    contour,reform(bbr[q,*,0,*]),pror1,vecz,nlev=20,color=col(q),/fill,title=string(q)
;;   wait,0.12
;; endfor
;
end

