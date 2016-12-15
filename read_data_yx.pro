pro read_data_yx, nm, yx, lm_tot
n = ' '
n_clu = 28
openr, 1, 'yx_vs_m_500.dat'
name = strarr(n_clu)  
yx = fltarr(n_clu)
e_z = fltarr(n_clu)
lm_tot = fltarr(n_clu)

for i=0,n_clu-1 do begin
    readf, 1, n, p1, p2, p3, format='(A19,1X,F7.6,6X,F7.5,6X,F7.5)'
    name[i] = n
    lm_tot[i] = p2
    yx[i]  = p3
    e_z[i] = p3
endfor
close, 1

nm = name

!psym = 4
plot,lm_tot,yx,xtitle='Log(Mtot)',ytitle='Yx',yrange=[13.0,16.0]
!psym=0
;oplot,lm_tot,yx_fit
;
;plot,yx,yx_fit
;stop

SIXLIN,lm_tot,yx,a,asig,b,bsig

;print,a,asig

;print,b,bsig

;y1 = b[1]*lm_tot + a[1]

;oplot,lm_tot,y1

; for each data point est error and sum chisq
for j = 1,20 do begin
chisq = 0.
yx_err = yx * .001*j ; 1% errors
for i = 1,n_clu - 1 do begin
ypred = b[2]*lm_tot[i] + a[2]
chisq = chisq + ((yx[i] - ypred)^2)/yx_err[i]^2

endfor

chisq = (chisq/(n_clu - 1))
print,'chisq:',chisq, '   Err:',0.01*j
endfor
;stop
end
