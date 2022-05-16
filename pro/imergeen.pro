;#; program to merge the dat/imf_{b,v}profilesq_?.sav files in one single sav containing the profiles for the B,v fields
pro imergeen,path,str

    path0 = path
    path = path

;#; understand which kind of energy is being computed'
case str of
 'b': mssg='Merging magnetic energy files'
 'v': mssg='Merging kynetic energy files'
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
        print,'i',i
        file = path+string(i,format=fileformat)
        restore,file
    
        if (i eq 1) then begin
         case str of
          'b': b_en = temporary(total_en)
          'v': v_en = temporary(total_en)
         endcase
         case str of
          'b': b_comp = temporary(comp) 
          'v': v_comp = temporary(comp)
         endcase
          tt = temporary(t)
            
        endif else begin

         case str of
          'b': b_en = [b_en,temporary(total_en)]
          'v': v_en = [v_en,temporary(total_en)]
         endcase
         case str of
          'b': b_comp = [b_comp,temporary(comp)]  
          'v': v_comp = [v_comp,temporary(comp)]
         endcase
          tt = [tt,temporary(t)]

        endelse
          
        t = 0
        total_en = 0
        comp = 0
        comp2 = 0
;#; end of the main for cycle
    endfor
         case str of
          'b': mssg = 'Sto salvando il file.. attendere.. b_energy.sav'
          'v': mssg = 'Sto salvando il file.. attendere.. v_energy.sav'
         endcase
         case str of
          'b': save_file = path+'dat/b_energy.sav'
          'v': save_file = path+'dat/v_energy.sav'
         endcase
         case str of
          'b': save,tt,b_comp,b_en,filename=save_file
          'v': save,tt,v_comp,v_en,filename=save_file
         endcase

   
end
