;#; program to merge the dat/imf_{b,v}profilesq_?.sav files in one single sav containing the profiles for the B,v fields
pro imergevmf,path,d3,mf

    path0 = path
    path = path

;    restore,path+'dat/spectrum.sav'

result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,d3,my,mm,nanf,nz,nzcon,i2d,m2d,n2d
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
    jmin = 0
    jmax = 1

    save_file = path+'dat/imf_vprofiles.sav'
; scopro il numero di sezioni sav
    fileformat1 = '("dat/imf_vprofiles_",I0,".sav")'
    nsectionsdat = find_nsectionsdat(path)
    nsectionssav = find_nsectionssav(path)
;#; Marco, check wheter imf_bprofiles.sav exists
    ex_0 = file_test(save_file)
    if (ex_0 eq 0) then begin
     goto, create_merge
    endif else begin
;#; Marco, aggiungere check sul fatto che imf sia già uptodate!
;     print,nsectionsdat,nsectionssav
;     print,save_file
;     print,path+string(nsectionssav,format=fileformat1)
;     cmpt = test_up_to_date2(path,string(nsectionssav,format=fileformat1))
     cmpt = test_up_to_date_files(save_file,path+string(nsectionssav,format=fileformat1))
     print,cmpt
     if (cmpt eq 0) then begin
      print,'no need to merge the files'
      goto, endprogram
     endif else begin
      goto, add_merge
     endelse
    endelse

    if (nsectionssav eq 1) then begin
;#; no need to merge, just mv the file
     spawn, "mv " +strtrim(string(nsectionssav,format=fileformat1))+" "+save_file
     goto, endprogram
    endif

    create_merge: print,'non esiste il file, lo creo'
    add_merge: print,'unisco i nuovi files'
    for i=1,nsectionsdat do begin
     dummm=fix(file_test(path+string(i,format=fileformat1)))
     if (dummm eq 0) then begin
      call_procedure, 'ireadsc2',path,nsectionssav,mf,d3
     endif
    endfor
    nsections = 1

    while (file_test(path+string(nsections,format=fileformat1))) do begin
        print,'guardo se esiste la sezione ',nsections,file_test(path+string(nsections,format=fileformat1))
        nsections = nsections + 1
    endwhile
    nsections = nsections - 1
    print,'nsections=',nsections
    if (nsections lt nsectionsdat) then begin
     print,'qualcosa non va, mancano i sav',path,nsections,nsectionsdat
    endif

    nt = ulonarr(nsections)
    fileformat = '("dat/imf_vprofiles_",I0,".sav")'
    for i = 1,nsections do begin
        print,'i',i
        file = path+string(i,format=fileformat)
        restore,file
    
        if (i eq 1) then begin
            vr = temporary(mf_r)
            vt = temporary(mf_t)
            vz = temporary(mf_z)
            tt = temporary(t)

            print,'memory1',memory(/current)
        endif else begin

            vz=[vz,temporary(mf_z)]
            vt=[vt,temporary(mf_t)]
            vr=[vr,temporary(mf_r)]
            tt = [tt,temporary(t)]
            print,'memory2',memory(/current)

        endelse

        mf_r = 0
        mf_t = 0
        mf_z = 0
        t = 0
;#; end of the main for cycle
    endfor
    lx = 100.d
    pror1 = findgen(lx+1)/lx*1.d
    pror2 = findgen(lx+2)/lx*1.d - 1.d/(2.*lx)
        
   print,'Sto salvando il file.. attendere.. imf_vprofiles.sav'


   save,tt,pror1,pror2,vr,vt,vz,filename=save_file
;    for i = 1,nsections do begin
;        print,'removing the partial files',i
;        file = path+string(i,format=fileformat)
;        spawn,'rm -f '+file
;    endfor
   
   endprogram: print,'sono alla fine, ho saltato unione dei files'
end
