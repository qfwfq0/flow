;#; this routines check wheter a directory exists. If it doesn't then it creates.

;#; written by Marco Veranda, 05.09.2014

pro check_dir,dir

 print,dir
;#;
 chk = file_test(dir,/directory)

 if (chk ne 1) then begin
  spawn,'mkdir -p '+dir
 endif

end
