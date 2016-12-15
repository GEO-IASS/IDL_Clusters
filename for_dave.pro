pro for_dave
restore, 'mass_temp_res_corr.sav'

omega_m = 0.3 & omega_lam = 0.7
e_z =  sqrt(omega_m*(1.0 + zed)^3 + omega_lam)
n = n_elements(zed)


for i=1,2 do begin
     if (i eq 1) then begin
;Do with r2500 - A2163
         x = where((mass_new[*,1,0,0] gt 0.0) and $
                   (strcompress(clu_name_mt, /remove_all) ne 'A2163'),nx)
         openw, 1, 'for_dave_2500_n.dat'
         printf, 1, 'R2500 measurements without A2163'
         r = 0
     endif

     if (i eq 2) then begin
;Do with r2500 - A2163
         x = where(mass_new[*,1,1,0] gt 0.0,nx)
         openw, 1, 'for_dave_500_n.dat'
         printf, 1, 'R500 measurements'
         r = 1
     endif

j=2

;Mass total
     mt = mass_new[x,1,r,0]
     yup = mass_new[x,1,r,2]
     ydown = mass_new[x,1,r,1]
     mterr = 0.5*abs(yup - ydown)

;Temp
     t = temps[x,j,r,0]
     xup = temps[x,j,r,2]
     xdown = temps[x,j,r,1]
     terr = 0.5*abs(xup - xdown)

;Mg
     mg = mass_new[x,0,r,0]
     yup = mass_new[x,0,r,2]
     ydown = mass_new[x,0,r,1]
     mgerr = 0.5*abs(yup - ydown)

;Y
     yx = alog10(temps[x,j,r,0]) + mass_new[x,0,r,0]
     xup = alog10(temps[x,j,r,2]) + mass_new[x,0,r,2]
     xdown = alog10(temps[x,j,r,1]) + mass_new[x,0,r,1]
     yxerr = 0.5*abs(xup - xdown)


     for k=0,nx-1 do printf, 1, clu_name_mt[x[k]], e_z[x[k]], t[k], terr[k], $
       yx[k], yxerr[k], mg[k], mgerr[k], mt[k], mterr[k]

     close, 1
endfor

end
