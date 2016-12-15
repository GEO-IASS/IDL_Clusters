pro read_data_temp, nm, parm
n = ' '
n_clu = 71
openr, 1, 'final_temp.dat'
name = strarr(n_clu) & parms = fltarr(n_clu,8)

for i=0,n_clu-1 do begin
    readf, 1, n, p1, p2, p3, p4, format='(A14,5X,F14,F15,3X,F14,2X,F12)'
    readf, 1, p5, p6, p7, p8
    readf, 1, f1, f2, f3
    name[i] = n
    parms[i,*] = [p1,p2,p3,p4,p5,p6,p7,p8]
endfor
close, 1

nm = name
parm = parms

end
