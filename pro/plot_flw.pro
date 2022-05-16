pro plot_flw, path, d3, itp,  hel

   RR = read_major_radius(path)
letter="143B;"
greekchi='!9' + String(letter) + '!X'
;path = './'
;d3 = 0
;hel = -10.

outfile1 = path+'/dat/itp/'+strtrim(string(itp,format='(i0)'))+'/flowxy.sav'
outfile2 = path+'/dat/itp/'+strtrim(string(itp,format='(i0)'))+'/chi'+strtrim(string(hel,format='(i0)'))+'.sav'
chk1 = file_test(outfile1)
chk2 = file_test(outfile2)

if (chk1 ne 1) then begin
 call_procedure,'floxwxy2',path,itp,ieli,d3,nth,nzz
endif
if (chk2 ne 1) then begin
 nth = 128
 nzz = 512
 call_procedure,'idl_savechi',path,itp,ieli,d3,nth,nzz
endif
restore,outfile1
restore,outfile2

xx=pror#cos(th)
yy=pror#sin(th)

zeta = 0

filename=path+'dat/itp/'+strtrim(string(itp,format='(i0)'))+'/plot_flowxy.eps'
phd_graphopener,filename=filename,xsize=10.,ysize=10.

contour,rchi[*,*,zeta],xx,yy,nlev=10,/iso,position=[0.15,0.15,0.92,0.92],xtitle='x',ytitle='y',title=greekchi+'  equilibrium and m=1 n=hel flow'
partvelvec2, vvxx,vvyy,vxxx,vyyy,length=0.04,color=250,veccolors=100,fraction=0.05,/over

phd_graphcloser,filename
end
