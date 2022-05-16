function cross_product,v1r,v1t,v1z,v2r,v2t,v2z
;#; this function computes the cross product of 2 arrays, v1,v2 
;v1, v2 are vectors in the real space


 array = make_array(n_elements(v1r),3,/double,value=0.d)

;#; component r 
 array[*,0] = v1t*v2z - v1z*v2t
 array[*,1] = v1z*v2r - v1r*v2z
 array[*,2] = v1r*v2t - v1t*v2r
 
 return,array

end
