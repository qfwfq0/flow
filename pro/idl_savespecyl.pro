pro idl_savespecyl, path

 print,path
;#; this program reads the output of SpeCyl and changes it creating an input file that can be used to restart SpeCyl
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

  print,'saving at the time',save_time,'  of section',nspc
  idl_file=path+'dat/specyl_'+strtrim(string(nspc,format='(i0)'))+'_all.sav'
  print,'restoring...',idl_file
  restore,idl_file
  
  help,br
  dimension = size(br)
  if (dimension[0] eq 3) then begin ;#; 2D sim
   br = reform(br[save_time,*,*])
   bt = reform(bt[save_time,*,*])
   bz = reform(bz[save_time,*,*])
   vr = reform(vr[save_time,*,*])
   vt = reform(vt[save_time,*,*])
   vz = reform(vz[save_time,*,*])
   
;#; here I insert a procedure that modifies the field!
;#; first of all: phase adding
;#; modes I want to modify
   min = 1
   max = 1
;#; phase to be added to the mode
   ph = 1.0d
   modif= 'r'
   brr=idl_add_phase(br,modif,ph,min,max)
   modif= 't'
   btt=idl_add_phase(bt,modif,ph,min,max)
   modif= 'z'
   bzz=idl_add_phase(bz,modif,ph,min,max)


   comment_file=path+'dat/specyl_'+strtrim(string(nspc,format='(i0)'))+'_end_com.dat'
   save_file = path+'dat/specyl_'+strtrim(string(nspc,format='(i0)'))+'_end_mod.dat'
   openw,unitc,comment_file,/get_lun
    printf,unitc,'Added a phase '+strtrim(string(ph,'(f13.5)'))+' rad to the modes between  '+strtrim(string(min,'(i0)'))+' and '+strtrim(string(min,'(i0)'))+' .'
   free_lun,unitc

   openw,unitw,save_file,/get_lun,/f77_unformatted
;    writeu, unitw, long64(itp), long64(n_elements(br[*,0])),long64(n_elements(br[0,*]))
    writeu, unitw, t(save_time), n_elements(br[*,0]),n_elements(br[0,*])
    
     for j = 0, long64(n_elements(br[0,*])) - 1 do begin
      writeu, unitw, vt[0,j], vz[0,j], bt[0,j], bz[0,j]
      for l = 0, long64(n_elements(br[*,0])) - 1 do begin
       writeu, unitw, vr[l,j], vt[l,j], vz[l,j], br[l,j], bt[j,j], bz[l,j]
      endfor
     endfor
   free_lun,unitw
    print, t(save_time), n_elements(br[*,0]),n_elements(br[0,*])
    help,bt
  endif else begin ;#; 3D SIM
  endelse

end
