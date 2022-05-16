pro idl_energy,path

print,path
nsectionssav = find_nsectionssav(path)
fileformat1 = '("dat/specyl_",I0,"_all.sav")'
fileformatb = '("dat/b_energy_",I0,".sav")'
fileformatv = '("dat/v_energy_",I0,".sav")'

savfile = 'dat/specyl_'+strtrim(string(nsectionssav))+'_all.sav'
en_file = 'dat/b_energy.sav'
uptodate =test_up_to_date_files(en_file,savfile)
en_test = file_test(en_file)

print,uptodate,en_test,nsectionssav
;stop
if (uptodate eq 1 or en_test eq 0) then begin
for j = 1, nsectionssav do begin

 print,'file: ',path+string(j,format=fileformat1)
 file=path+string(j,format=fileformat1)
 bfile = path+string(j,format=fileformatb)
   ex_file = file_test(bfile)   
   if (ex_file ne 1) then begin
    print,'need to create the energy files' 
    restore,file
;#; compute magnetic energy
    call_procedure,'energy',t,br,bt,bz,'b',j,path
;#; compute kynetic energy
    call_procedure,'energy',t,vr,vt,vz,'v',j,path
   endif else begin
    print,'already exists'
    continue
   endelse
endfor

 call_procedure,'imergeen',path,'b'
 call_procedure,'imergeen',path,'v'
endif
 call_procedure,'idelen',path,'b'
 call_procedure,'idelen',path,'v'

end
