pro idl_savepixie, path

;#; this program reads the output of SpeCyl and saves it in a format that will be read by PIXIE3D;
;#; REMARK: only the perturbations to the equilibrium field are written.

 nsectionssav = find_nsectionssav(path)
 xxx = read_settings(path,nsectionssav)

 if (nsectionssav ne n_elements(xxx.computed)) then begin
  print,'produce the correct specyl_?_all.sav files!'
 endif

;#; I want to dump on disk the time xxx.max
  itp = xxx.max
;#; this itp is written on the file specyl_?_all.sav with ?=
  nspc = min(where( itp gt xxx.cumulated )) + 1 + 1
      print,nspc
      print,xxx.cumulated
;#;at the time 
  save_time = itp - xxx.cumulated(nspc - 1 - 1) -1 

  print,'saving at the time',save_time,'of section',nspc
  idl_file=path+'dat/specyl_'+strtrim(string(nspc,format='(i0)'))+'_all.sav'
  print,'restoring...',idl_file
  restore,idl_file
  
  help,br
  dimension = size(br)
  if (dimension[0] eq 3) then begin ;#; 2D sim
   br = reform(br[save_time,*,1:*])
   bt = reform(bt[save_time,*,1:*])
   bz = reform(bz[save_time,*,1:*])
;   vr = reform(vr[save_time,*,*])
;   vt = reform(vt[save_time,*,*])
;   vz = reform(vz[save_time,*,*])
;#; now I interpol the field on a radial-pixie-like-mesh made of 128+2 elements.
;   nxd = long64(128)
;   r1 = findgen(dimension(2))/(dimension(2)-1)
;   r2 = findgen(dimension(2)+1)/(dimension(2)-1) - 1./(2.*(dimension(2)-1))
;   p3d_r = findgen(nxd+2)/nxd - 1./(2.d*nxd)
;
;   rbr2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   ibr2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   rbt2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   ibt2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   rbz2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   ibz2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   rvr2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   ivr2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   rvt2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   ivt2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   rvz2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   ivz2 = make_array(nxd+2,dimension(3),/double,value=0.)
;   for q = 0, dimension(3) - 1 do begin
;    rbr2(*,q) = interpol(real_part(br[*,q]),r1,p3d_r,/spline)
;    ibr2(*,q) = interpol(imaginary(br[*,q]),r1,p3d_r,/spline)
;    rbt2(*,q) = interpol(real_part(bt[*,q]),r2,p3d_r,/spline)
;    ibt2(*,q) = interpol(imaginary(bt[*,q]),r2,p3d_r,/spline)
;    rbz2(*,q) = interpol(real_part(bz[*,q]),r2,p3d_r,/spline)
;    ibz2(*,q) = interpol(imaginary(bz[*,q]),r2,p3d_r,/spline)
;    rvr2(*,q) = interpol(real_part(vr[*,q]),r1,p3d_r,/spline)
;    ivr2(*,q) = interpol(imaginary(vr[*,q]),r1,p3d_r,/spline)
;    rvt2(*,q) = interpol(real_part(vt[*,q]),r2,p3d_r,/spline)
;    ivt2(*,q) = interpol(imaginary(vt[*,q]),r2,p3d_r,/spline)
;    rvz2(*,q) = interpol(real_part(vz[*,q]),r2,p3d_r,/spline)
;    ivz2(*,q) = interpol(imaginary(vz[*,q]),r2,p3d_r,/spline)
;   endfor
;   br2 = dcomplex(rbr2,ibr2)
;   bt2 = dcomplex(rbt2,ibt2)
;   bz2 = dcomplex(rbz2,ibz2)
;   vr2 = dcomplex(rvr2,ivr2)
;   vt2 = dcomplex(rvt2,ivt2)
;   vz2 = dcomplex(rvz2,ivz2)
   
   save_file = path+'dat/spc_fft_t'+strtrim(string(itp,format='(i0)'))+'.bin'
   openw,unitw,save_file,/get_lun,/f77_unformatted
    writeu,unitw, long64(n_elements(br[*,0]))
    writeu,unitw, long64(n_elements(br[0,*]))
    writeu,unitw,long64(1)
    writeu,unitw,br
    writeu,unitw,bt
    writeu,unitw,bz
;    writeu,unitw,vr2
;    writeu,unitw,vt2
;    writeu,unitw,vz2
   free_lun,unitw
    print, long64(n_elements(br[*,0]))
    print, long64(n_elements(br[0,*]))
    help,bt
  endif else begin ;#; 3D SIM
  endelse

end
