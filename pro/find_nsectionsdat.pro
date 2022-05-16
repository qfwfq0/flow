function find_nsectionsdat,path
    fileformat1 = '("dat/specyl_",I0,"_all.dat")'
    nsectionsdat = 1
;        print,'guardo se esiste la sezione ',nsectionsdat
;        print,path+string(nsectionsdat,format=fileformat1)
    while (file_test(path+string(nsectionsdat,format=fileformat1))) do begin
        print,'guardo se esiste la sezione *dat',nsectionsdat
;        print,'file: ',path+string(nsectionsdat,format=fileformat1)
;        print,'risultato: ',file_test(path+string(nsectionsdat,format=fileformat1))
        nsectionsdat = nsectionsdat + 1
    endwhile
;    print,'nsectionsdat',nsectionsdat
    nsectionsdat = nsectionsdat - 1
    return,nsectionsdat
end
