function 2dconvo,mcap,ncap
;#; function that computes all the janz of the modes that make convolution with the mode mcap,ncap
;#; 2D case, very easy:
;#; a) the mode (mcap,ncap) makes convolution with itself
;#; b) the mode (0,0) makes convolution with itself

path='./'
;input='spectrum.sav'
;restore,path+'dat/'+strtrim(input)
result = file_test(path+'dat/spectrum.sav')
if (result ne 1) then begin
 call_procedure, 'read_spectrum',path,my,mm,nanf,nz
endif else begin
 restore,path+'dat/spectrum.sav'
endelse

mcap=1
ncap=-10
;#; consistency check.
if (mcap gt my) then begin
 print,'wrong m',mcap
 stop
endif
if (ncap lt nanf(mcap) or ncap gt nanf(mcap)+nz(mcap)-1) then begin
 print,'wrong n',ncap
 stop
endif

janzf=[0,janz(mcap,ncap,my,mm,nanf,nz)]

return,janzf
end
