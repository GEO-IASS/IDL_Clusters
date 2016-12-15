pro read_masspar_dat, nm, parm
n = ' '
n_clu = 71
close,1
openr, 1, 'xmm_masspar.dat'
name = strarr(n_clu) & parms = fltarr(n_clu,8)

for i=0,n_clu-1 do begin
    readf, 1, n, p1, p2, p3, p4, p5, p6, p7
    name[i] = n
    parms[i,*] = [n, p1,p2,p3,p4,p5,p6,p7]
endfor
close, 1

nm = name
parm = parms

end
