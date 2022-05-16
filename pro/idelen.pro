;#; program to merge the dat/imf_{b,v}profilesq_?.sav files in one single sav containing the profiles for the B,v fields
pro idelen,path,str

    path0 = path
    path = path

;#; understand which kind of energy is being computed'
case str of
 'b': mssg='Deleting magnetic energy files'
 'v': mssg='Deleting kynetic energy files'
endcase
; scopro il numero di sezioni sav
    fileformat1 = '("dat/'+str+'_energy_",I0,".sav")'
    nsections = 1
    while (file_test(path+string(nsections,format=fileformat1))) do begin
;        print,'guardo se esiste la parte ',nsections
        nsections = nsections + 1
    endwhile
    nsections = nsections - 1
;    print,'nsections=',nsections
    nsectionsdat = find_nsectionsdat(path)
    if (nsections lt nsectionsdat) then begin
     print,'qualcosa non va, mancano i sav'
     stop
    endif

    nt = ulonarr(nsections)
    fileformat = '("dat/'+str+'_energy_",I0,".sav")'
    for i = 1,nsections do begin
        file = path+string(i,format=fileformat)
        print,'i',i,file
        spawn,'rm -f '+file
    endfor

   
end
