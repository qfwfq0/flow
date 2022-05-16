function nonlinearconvo,path,vect1,vect2,j,d3
;#; this function computes the modes resulting from the nonlinear coupling of two quantities, function of the radius, vect1[radius,modes] and vect2[radius,modes]

;path='./'
restore,path+'dat/spectrum.sav'
 array = make_array(n_elements(vect1[*,0]),/dcomplex,value=0.d)
     mn = mnum(j,nz,nanf,d3)
     cnvl = convo(path,mn[0],mn[1],d3)
     help,cnvl
    for p = 0, n_elements(cnvl[0,*])-1 do begin
     j1 = cnvl[0,p]
     j2 = cnvl[1,p]
     mn1 = mnum(j1,nz,nanf,d3)
     mn2 = mnum(j2,nz,nanf,d3)

     print,'nonlinearconvo0',p,j,j1,j2
;#; definition of the multiplication factor
     fact = mlpf(mn1,mn2)

     array = array + fact*vect1[*,j1]*vect2[*,j2] 
;#; Marco, end of the convolution cycle
    endfor
    return,array
end
