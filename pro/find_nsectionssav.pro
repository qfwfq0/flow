function find_nsectionssav,path
    fileformat1 = '("dat/specyl_",I0,"_all.sav")'
    nsectionssav = 1
    while (file_test(path+string(nsectionssav,format=fileformat1))) do begin
;        print,'guardo se esiste la sezione ',nsectionssav
;        print,'file: ',path+string(nsectionssav,format=fileformat1)
;        print,'risultato: ',file_test(path+string(nsectionssav,format=fileformat1))
        nsectionssav = nsectionssav + 1
    endwhile
;    print,'nsectionssav',nsectionssav
    nsectionssav = nsectionssav - 1
    return,nsectionssav
end
