function integral2,vect,xx
;#; general purpose function to calculate the integral of a function (vect) written on the mesh x2 (vect)
;#; x2=findgen(lx+2)/lx-1./(2.*lx)
;#; the result will be stored in the mesh x1(i) = findgen(lx+1)/lx

 integ = make_array(n_elements(xx)-1,n_elements(vect(0,*)),/double,value=0.)
 yy = findgen(n_elements(xx)-1) / (n_elements(xx)-2)
 dx = 1.d / (n_elements(xx)-2)
 
 for q=0, n_elements(vect(0,*))-1 do begin
  integ(0,q) = 0.d
  for i=1,n_elements(yy)-1 do begin
;   integ(i,q) = integ(i-1,q) + vect(i,q) * dx
   integ(i,q) = integ(i-1,q) + vect(i-1,q) * dx
  endfor
 endfor
 return,integ
end

function integralc,vect,xx
;#; general purpose function to calculate the integral of a function (vect) along the mesh x2 (vect)
;#; x2=findgen(lx+2)/lx-1./(2.*lx)
;#; the result will be stored in the mesh x1(i) = findgen(lx+1)/lx

 integ = make_array(n_elements(xx)-1,/dcomplex,value=0.)
 yy = findgen(n_elements(xx)-1) / (n_elements(xx)-2)
 dx = 1.d / (n_elements(xx)-2)
 
; for q=0, n_elements(vect(*))-1 do begin
  integ(0) = 0.d
  for i=1,n_elements(yy)-1 do begin
   integ(i) = integ(i-1) + vect(i-1) * dx
  endfor
; endfor
 return,integ
end

pro idl_savea,path,it,d3
;#; this routine computes the vector potential components useful to the helical flux function.
;#; it computes Az=-integ(dr*Btheta) and r*Atheta=integ(dr*r*Bz)
;#; Marco, 06 settembre 2019, modifico per calcolare A in tutti i tempi (nel caso it=-2)

;    path = path+'/dat'
;  path='./dat/'
  RR = read_major_radius(path)

;  restore,path+'dat/spectrum.sav'
  result = file_test(path+'dat/spectrum.sav')
  if (result ne 1) then begin
   call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
  endif else begin
   restore,path+'dat/spectrum.sav'
   nzcon=nz
   nzcon[0]=nz[0]+1
  endelse

  nsectionsdat = find_nsectionsdat(path)
  print,nsectionsdat
  xxx = read_settings(path,nsectionsdat)
 
  print,it,' 3d=',d3
  help,xxx,/str
  idl_file=path+'dat/ibprofiles.sav'
;#; attenzione in ibprofiles.sav c'è già il fattore 2 che serve per fare una corretta antifft
  dumdum = file_test(idl_file)
  if (dumdum ne 1) then begin
   print,'ibprofiles not present, making it'
   call_procedure, 'isavemodeprofiles',path,0,d3
  endif else begin
;#; check whether it's up to date
   nsectionssav = find_nsectionssav(path)
   savfile = 'dat/specyl_'+strtrim(string(nsectionssav))+'_all.sav'
   uptodate =test_up_to_date_files(idl_file,savfile)
   if (uptodate eq 1) then begin
    print,'ibprofiles to be updated, making it'
    call_procedure, 'isavemodeprofiles',path,0,d3
   endif
  endelse
  print,'restoring...',idl_file
  restore,idl_file

  if (it eq -2) then begin
 ;#; Marco compute the FOurier components of vector potential at all times
 ;#; I declare the vector potential real and imaginary part.
 ;#; REMIND: guage assumption, Ar=0.
;   rpdaz = reform(real_part(pdbz[*,*,*])) * 0.d
;   rpdat = reform(real_part(pdbz[*,*,*])) * 0.d
;   ipdaz = reform(real_part(pdbz[*,*,*])) * 0.d
;   ipdat = reform(real_part(pdbz[*,*,*])) * 0.d
   at = make_array(n_elements(pdbz[*,0,0]),n_elements(pdbz[0,*,0])-1,n_elements(pdbz[0,0,*]),/dcomplex,value=0.)
   az = make_array(n_elements(pdbz[*,0,0]),n_elements(pdbz[0,*,0])-1,n_elements(pdbz[0,0,*]),/dcomplex,value=0.)
   pdat = make_array(n_elements(pdbz[*,0,0]),n_elements(pdbz[0,*,0]),n_elements(pdbz[0,0,*]),/dcomplex,value=0.)
   pdaz = make_array(n_elements(pdbz[*,0,0]),n_elements(pdbz[0,*,0]),n_elements(pdbz[0,0,*]),/dcomplex,value=0.)
   uno = make_array(n_elements(reform(real_part(pdbz[0.,*,*]))),/double,value=1.)
   pror3 = pror2#uno
 
   for dit=0,n_elements(pdbz[*,0,0])-1 do begin
;#; ciclo sui modi
    for q=0, n_elements(pdbz[dit,0,*]) - 1 do begin
     if (q eq 0) then begin
;#; Marco, 05 dic 2019: avevo calcolato l'integrale per il calcolo di A_theta ma dimenticato di dividere per il raggio
      at[dit,*,q] = integralc( reform( pdbz[dit,*,q] ) * pror2,pror2 ) / pror2
      az[dit,*,q] = integralc( - reform( pdbt[dit,*,q] ) ,pror2 )
     endif else begin
 ; #; divido per due perché in ibprofiles il modo è già stato moltiplicato per 2
      at[dit,*,q] = integralc( reform( pdbz[dit,*,q] / 2.d ) * pror2,pror2 ) / pror2
      az[dit,*,q] = integralc( - reform( pdbt[dit,*,q] / 2.d) ,pror2 )
     endelse
    endfor
   endfor

;#; interpolate in mesh X2
   for it=0,n_elements(at[*,0,0])-1 do begin
    for j=0,n_elements(at[0,0,*])-1 do begin
    mn = mnum(j,nzcon,nanf,d3)
    pdat[it,*,j] = interpol(reform(at[it,*,j]), pror1, pror2 )
    pdaz[it,*,j] = interpol(reform(az[it,*,j]), pror1, pror2 )
    if (mn[0] eq 0) then begin
     pdat[it,0,j] = -pdat[it,1,j]
     pdaz[it,0,j] = -pdaz[it,1,j]
    endif else begin
    if (mn[1] eq 1) then begin
     pdat[it,0,j] = pdat[it,1,j]
     pdaz[it,0,j] = pdaz[it,1,j]
    endif else begin
     pdat[it,0,j] = -pdat[it,1,j]
     pdaz[it,0,j] = -pdaz[it,1,j]
    endelse
    endelse
    endfor
   endfor
   save_file = path+'dat/iaprofiles.sav'
   tt = pdt
   help,pdat
   save,pror1,pror2,tt,pdaz,pdat,filename=save_file

  endif else begin
   if (it lt xxx.min) then begin
    print,'The time requested is not present on this ibprofiles.sav. Remake it with the correct time-slice selection!'
    stop
   endif
   
   dit = fix(it - xxx.min)
   help,pdbz
   print,'it0= ',it,' dit0= ',dit 
   if (dit ge n_elements(pdbz[*,0,0])) then begin
    dit = n_elements(pdbz[*,0,0])-1
   endif
   print,'it= ',it,' dit= ',dit 
 ;#; I declare the vector potential real and imaginary part.
 ;#; REMIND: guage assumption, Ar=0.
   rpdaz = reform(real_part(pdbz[dit,*,*])) * 0.d
   rpdat = reform(real_part(pdbz[dit,*,*])) * 0.d
   ipdaz = reform(real_part(pdbz[dit,*,*])) * 0.d
   ipdat = reform(real_part(pdbz[dit,*,*])) * 0.d
   dat = make_array(n_elements(pdbz[dit,*,0])-1,n_elements(pdbz[dit,0,*]),/dcomplex,value=0.)
   daz = make_array(n_elements(pdbz[dit,*,0])-1,n_elements(pdbz[dit,0,*]),/dcomplex,value=0.)
 help,daz
   uno = make_array(n_elements(reform(real_part(pdbz[dit,*,*]))),/double,value=1.)
   pror3 = pror2#uno
 
 
   print,'integrating real(Az), itp='+it+','
   rpdaz = integral2(-reform(real_part(pdbt[dit,*,*])),pror2)
   print,'integrating im(Az)'
   ipdaz = integral2(-reform(imaginary(pdbt[dit,*,*])),pror2)
;#; Marco, 05 dic 2019: avevo calcolato l'integrale per il calcolo di A_theta ma dimenticato di dividere per il raggio
   print,'integrating real(r*At)'
   rpdat = integral2(reform(real_part(pdbz[dit,*,*])) * pror3,pror2)/pror2
   print,'integrating im(r*At)'
   ipdat = integral2(reform(imaginary(pdbz[dit,*,*])) * pror3,pror2)/pror2
   for q=0, n_elements(pdbz[dit,0,*]) - 1 do begin
    if (q eq 0) then begin
     dat[*,q] = integralc( reform( pdbz[dit,*,q] ) * pror2,pror2 )/pror2
     daz[*,q] = integralc( - reform( pdbt[dit,*,q] ) ,pror2 )
    endif else begin
 ;#; divido per due perché in ibprofiles il modo è già stato moltiplicato per 2
     dat[*,q] = integralc( reform( pdbz[dit,*,q] / 2.d ) * pror2,pror2 )/pror2
     daz[*,q] = integralc( - reform( pdbt[dit,*,q] / 2.d) ,pror2 )
    endelse
   endfor
   print,'it2= ',it
   save_file = path+'dat/itp/'+strtrim(string(it,format='(i0)'))+'/acyl.sav'
   save,pror1,pror2,it,rpdaz,ipdaz,rpdat,ipdat,daz,dat,filename=save_file

  endelse
end
