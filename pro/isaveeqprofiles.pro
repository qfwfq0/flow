;#; program to merge the dat/bv_eq_?.sav files in one single sav containing the equilibrium profiles for the B,v fields
pro isaveeqprofiles,path

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
    fileformat1 = '("dat/bv_eq_",I0,".sav")'
    nsections = 1
    while (file_test(path+string(nsections,format=fileformat1))) do begin
        print,'guardo se esiste la sezione ',nsections
        nsections = nsections + 1
    endwhile
    nsections = nsections - 1
    print,'nsections=',nsections

    nt = ulonarr(nsections)
    fileformat = '("dat/bv_eq_",I0,".sav")'
    for i = 1,nsections do begin
        print,'i',i
        file = path+string(i,format=fileformat)
        restore,file
    
        if (i eq 1) then begin
            eqbz = temporary(bz0)
            eqbt = temporary(bt0)
            eqbr = temporary(br0)
            eqvz = temporary(vz0)
            eqvt = temporary(vt0)
            eqvr = temporary(vr0)
            eqt = temporary(t)

            print,'memory1',memory(/current)
        endif else begin

            eqbz=[eqbz,temporary(bz0)]
            eqbt=[eqbt,temporary(bt0)]
            eqbr=[eqbr,temporary(br0)]
            eqvz=[eqvz,temporary(vz0)]
            eqvt=[eqvt,temporary(vt0)]
            eqvr=[eqvr,temporary(vr0)]
            eqt = [eqt,temporary(t)]
            print,'memory2',memory(/current)

        endelse

        bz0 = 0
        bt0 = 0
        br0 = 0
        vz0 = 0
        vt0 = 0
        vr0 = 0
        pt0 = 0
;#; end of the main for cycle
    endfor
    pror1 = findgen(lx+1)/lx*1.d
    pror2 = findgen(lx+2)/lx*1.d - 1.d/(2.*lx)
        
   print,'Sto salvando il file.. attendere.. ieqprofiles.sav'

   save_file = path+'dat/bv_eq.sav'
   save,eqt,pror1,pror2,eqbr,eqbt,eqbz,eqvz,eqvt,eqvr,filename=save_file
   
end
