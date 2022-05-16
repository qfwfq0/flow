function integral3,vect,xx
;#; general purpose function to calculate the integral of a function (vect) along the mesh x2 (vect)
;#; x2=findgen(lx+2)/lx-1./(2.*lx)
;#; the result will be stored in the mesh x1(i) = findgen(lx+1)/lx

;#; integral from the edge to the core.
 integ = make_array(n_elements(xx)-1,n_elements(vect(0,*)),/double,value=0.)
 yy = findgen(n_elements(xx)-1) / (n_elements(xx)-2)
 dx = 1.d / (n_elements(xx)-2)
 
 for q=0, n_elements(vect(0,*))-1 do begin
  integ(0,q) = 0.d
;  for i=1,n_elements(yy)-1 do begin
  for i=n_elements(yy)-2,0,-1 do begin
   integ(i,q) = integ(i+1,q) + vect(i+1,q) * dx
  endfor
 endfor
 return,integ
end
