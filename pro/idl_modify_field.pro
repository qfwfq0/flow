pro idl_modify_field,br,bt,bz,vr,vt,vz

  
  help,br
  dimension = size(br)
  if (dimension[0] eq 3) then begin ;#; 2D sim 
;#; here I insert a function that modifies the field!
;#; Marco, add a phase to the field   
  endif else begin ;#; 3D sim
  endelse

end
