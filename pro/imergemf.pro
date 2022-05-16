;#; program to merge the dat/imf_{b,v}profilesq_?.sav files in one single sav containing the profiles for the B,v fields
pro imergemf,path

    path0 = path
    path = path

;    restore,path+'dat/spectrum.sav'
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse
    jmin = 0
    jmax = 1

; scopro il numero di sezioni sav
    fileformat1 = '("dat/imf_vprofiles_",I0,".sav")'
    nsections = 1
    while (file_test(path+string(nsections,format=fileformat1))) do begin
        print,'guardo se esiste la sezione ',nsections
        nsections = nsections + 1
    endwhile
    nsections = nsections - 1
    print,'nsections=',nsections

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
            t = temporary(t)

            print,'memory1',memory(/current)
        endif else begin

            vz=[vz,temporary(mf_z)]
            vt=[vt,temporary(mf_t)]
            vr=[vr,temporary(mf_r)]
            t = [t,temporary(t)]
            print,'memory2',memory(/current)

        endelse

        mf_r = 0
        mf_t = 0
        mf_z = 0
        t = 0
;#; end of the main for cycle
    endfor
    pror1 = findgen(lx+1)/lx*1.d
    pror2 = findgen(lx+2)/lx*1.d - 1.d/(2.*lx)
        
   print,'Sto salvando il file.. attendere.. imf_vprofiles.sav'

   save_file = path+'dat/imf_vprofiles.sav'
   save,t,pror1,pror2,vr,vt,vz,filename=save_file
   
end
