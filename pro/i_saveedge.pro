;#; program to save the magnetic field at the edge
pro i_saveedge,path,d3,mf

; scopro il numero di sezioni sav
    nsectionsdat = find_nsectionsdat(path)
    nsectionssav = find_nsectionssav(path)
    fileformat1 = '("dat/imf_bprofiles_",I0,".sav")'
    edge_file = path+'dat/bfield_edge.sav'
    core_file = path+'dat/bfield_core.sav'
    ex_0 = file_test(edge_file)
    print,'aaa',ex_0
    if (ex_0 eq 0) then begin
     goto, create
    endif else begin
;     cmpt = test_up_to_date_files(edge_file,path+string(nsectionssav,format=fileformat1))
     cmpt = test_up_to_date_files(edge_file,path+'dat/imf_bprofiles.sav')
     if (cmpt eq 0) then begin
      print,'no need to merge the files'
      goto, endprogram
     endif else begin
      goto, create
     endelse
    endelse

    create: print,'creo edge-core files'
    save_file = path+'dat/imf_bprofiles.sav'
;#; Marco, check wheter imf_bprofiles.sav exists
     print,save_file
    ex_0 = file_test(save_file)
    if (ex_0 eq 0) then begin
     call_procedure,'imergebmf',path,d3,mf
    endif else begin
;#; Marco, aggiungere check sul fatto che imf sia già uptodate!
     cmpt = test_up_to_date_files(save_file,path+string(nsectionssav,format=fileformat1))
     print,cmpt
     if (cmpt eq 0) then begin
;      print,'no need to merge the files'
     endif else begin
      call_procedure,'imergebmf',path,d3,mf
     endelse
    endelse

    restore,save_file

    edgebr = reform(br[*,95:*,*,*])
    edgebt = reform(bt[*,95:*,*,*])
    edgebz = reform(bz[*,95:*,*,*])
    corebr = reform(br[*,0:6,*,*])
    corebt = reform(bt[*,0:6,*,*])
    corebz = reform(bz[*,0:6,*,*])
    pror1 = findgen(n_elements(br[0,*,0,0])) / (n_elements(br[0,*,0,0]) - 1)
    pror2 = findgen(n_elements(bt[0,*,0,0])) / (n_elements(br[0,*,0,0]) - 1) - 0.5 / (n_elements(br[0,*,0,0]) - 1)

   save,tt,pror1,pror2,edgebr,edgebt,edgebz,filename=edge_file
   save,tt,pror1,pror2,corebr,corebt,corebz,filename=core_file
   
   endprogram: print, 'tutto ok'
end
