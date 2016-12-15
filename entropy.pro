;;calc entropy at 500, 2500, calc K_500,adi
;;K_500,adi = 342 keVcm^2 (M500/10^14 msol)^(2/3.) E(z)^(-2/3)
;;h73^(-4/3)
pro entropy_dist_fig
omega_m = 0.3
omega_lam = 0.7
h0 = 70.0 ;km/s/Mpc

files ='a'
readcol, 'clu_500', files, format='A'
n_clu = n_elements(files)

set_plot, 'ps'
device, filename='entdist_fig.ps'
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0 
!p.charthick=2.0


;For first cluster
restore, strcompress(files[0],/remove_all)+'_dists.sav'
restore, strcompress(files[0],/remove_all)+'_derive.sav'
restore, strcompress(files[0],/remove_all)+'_basic.sav'

ez = sqrt(omega_m*(1.0 + z)^3 + omega_lam)
k500adi = 342.0 * (10.^(mass_tot_500-14.0))^(2./3.) * ez^(-2.0/3.) $
  * (h0/73)^(-4./3.)

plot, r_mpc/r500, entdist[*,0]/k500adi, /xlog, /ylog, xrange=[0.01,3], $
  xstyle=1, yrange=[0.02,5], ystyle=1, xtitle='Radius (r!B500!N)', $
  ytitle='K(r) / K!B500,adi!N', xmargin=[8,3]

for i=1,n_clu-1 do begin
    restore, strcompress(files[i],/remove_all)+'_dists.sav'
    restore, strcompress(files[i],/remove_all)+'_derive.sav'
    restore, strcompress(files[i],/remove_all)+'_basic.sav'

    ez = sqrt(omega_m*(1.0 + z)^3 + omega_lam)
    k500adi = 342.0 * (10.^(mass_tot_500-14.0))^(2./3.) * ez^(-2.0/3.) $
      * (h0/73)^(-4./3.)

    oplot, r_mpc/r500, entdist[*,0]/k500adi
endfor

rguide = (findgen(101)+0.005)*3.0/100
kguide = 1.40 * rguide^1.1

oplot, rguide, kguide, thick=3.0

device, /close
set_plot,'x'
end


;Calculate entropy values
pro get_entropy_vals

;2500
files = 'a'
readcol, '../clu_2500', files, format='A'
n_clu = n_elements(files)

for i=0,n_clu-1 do begin
    restore, '../'+strcompress(files[i],/remove_all)+'_basic.sav'
    restore, strcompress(files[i],/remove_all)+'_2500derive.sav_part'

    if ((temp_sample[0,0] eq -999) or (ne_sample[0,0] eq -999)) then begin
        print, 'Error analysis not possible, missing error information.' 
        radii_err = -999
    endif else begin
        n_temps = n_elements(temp_sample[*,0])
        n_nes = n_elements(ne_sample[*,0])

        res = [r2500, r2500_lim[0], r2500_lim[1]]

        ent_2500 = entropy_dist(r2500, ne_best_fit, temp_best_fit)
        entb = entropy_dist(res, ne_best_fit, temp_best_fit)

        ent_err = fltarr(3,n_nes,n_temps)

        for j=0,2 do begin
            for k=0,n_nes-1 do begin
                ne_pars = ne_sample[k,*]
                for l=0,n_temps-1 do begin
                    temp_pars = temp_sample[l,*]

                    ent_err[j,k,l] = entropy_dist(res[j], ne_sample[k,*], $
                                                     temp_sample[l,*])
                endfor
            endfor
        endfor

        ent_errtemp = fltarr(3,2)
        for j=0,2 do ent_errtemp[j,*] = get_limit(entb[j], ent_err[j,*,*], $
                                          0.8413)
        ent_2500_lim = fltarr(2)

        x = where(ent_errtemp[*,0] gt 0)
        ent_2500_lim[0] = min(ent_errtemp[x,0])
        ent_2500_lim[1] = max(ent_errtemp[*,1])

        save, name, ent_2500, ent_2500_lim, ent_err, $
          filename=strcompress(name,/remove_all)+'_2500ent.sav_part'
    endelse
endfor


;500
files = 'a'
readcol, '../clu_500', files, format='A'
n_clu = n_elements(files)

for i=0,n_clu-1 do begin
    restore, '../'+strcompress(files[i],/remove_all)+'_basic.sav'
    restore, strcompress(files[i],/remove_all)+'_500derive.sav_part'

    if ((temp_sample[0,0] eq -999) or (ne_sample[0,0] eq -999)) then begin
        print, 'Error analysis not possible, missing error information.' 
        radii_err = -999
    endif else begin
        n_temps = n_elements(temp_sample[*,0])
        n_nes = n_elements(ne_sample[*,0])

        res = [r500, r500_lim[0], r500_lim[1]]

        ent_500 = entropy_dist(r500, ne_best_fit, temp_best_fit)
        entb = entropy_dist(res, ne_best_fit, temp_best_fit)

        ent_err = fltarr(3,n_nes,n_temps)

        for j=0,2 do begin
            for k=0,n_nes-1 do begin
                ne_pars = ne_sample[k,*]
                for l=0,n_temps-1 do begin
                    temp_pars = temp_sample[l,*]

                    ent_err[j,k,l] = entropy_dist(res[j], ne_sample[k,*], $
                                                     temp_sample[l,*])
                endfor
            endfor
        endfor

        ent_errtemp = fltarr(3,2)
        for j=0,2 do ent_errtemp[j,*] = get_limit(entb[j], ent_err[j,*,*], $
                                          0.8413)
        ent_500_lim = fltarr(2)

        x = where(ent_errtemp[*,0] gt 0)
        ent_500_lim[0] = min(ent_errtemp[x,0])
        ent_500_lim[1] = max(ent_errtemp[*,1])

        save, name, ent_500, ent_500_lim, ent_err, $
          filename=strcompress(name,/remove_all)+'_500ent.sav_part'
    endelse
endfor

end


;K(r)/K_500,adi = 1.40(r/r500)^1.1
;do E(z)^(4/3) K (keVcm^2) vs kT for 2500, 500
;fit with BCES
;baseline

;Collect entropy values for 500 & 2500 to fit vs temp

pro collect_ent
omega_m = 0.3 & omega_lam = 0.7

;2500
files = 'a'
readcol, '../clu_2500', files, format='A'
n_clu = n_elements(files)
ey1ey2 = 0 

clu_names = strarr(n_clu)
ent_val = fltarr(3,n_clu)
temp_val = fltarr(3,n_clu)
zed = fltarr(n_clu)

openw, 1, 'ent_temp_2500.dat'
for i=0,n_clu-1 do begin
    restore, '../'+strcompress(files[i],/remove_all)+'_basic.sav'
    restore, strcompress(files[i],/remove_all)+'_2500derive.sav_part'
    restore, strcompress(files[i],/remove_all)+'_2500ent.sav_part'

    clu_names[i] = name
    ent_val[*,i] = [ent_2500, ent_2500_lim[0], ent_2500_lim[1]]
    temp_val[*,i] = [tx_2500, tx_2500_lim[0], tx_2500_lim[1]]
    zed[i] = z
    
;Need to print alog10(E^4/3 * K), err in that, alog10(T), err in that, dxdy

    e_z =  sqrt(omega_m*(1.0 + z)^3 + omega_lam)
    yaxis = alog10(e_z^(4./3.) * ent_2500)
    xaxis = alog10(tx_2500/5.0)
    xup = alog10(tx_2500_lim[1]/5.0) & xdown = alog10(tx_2500_lim[0]/5.0)
    yup = alog10(e_z^(4./3.) * ent_2500_lim[1]) 
    ydown = alog10(e_z^(4./3.) * ent_2500_lim[0]) 
    xerr = 0.5*abs(xup - xdown) & yerr = 0.5*abs(yup - ydown)

    printf, 1, xaxis, xerr, yaxis, yerr, ey1ey2

endfor
close, 1
save, clu_names, ent_val, temp_val, zed, filename='ent_2500.sav_part'

;500
files = 'a'
readcol, '../clu_500', files, format='A'
n_clu = n_elements(files)

clu_names = strarr(n_clu)
ent_val = fltarr(3,n_clu)
temp_val = fltarr(3,n_clu)
zed = fltarr(n_clu)

openw, 1, 'ent_temp_500.dat'
for i=0,n_clu-1 do begin
    restore, '../'+strcompress(files[i],/remove_all)+'_basic.sav'
    restore, strcompress(files[i],/remove_all)+'_500derive.sav_part'
    restore, strcompress(files[i],/remove_all)+'_500ent.sav_part'

    clu_names[i] = name
    ent_val[*,i] = [ent_500, ent_500_lim[0], ent_500_lim[1]]
    temp_val[*,i] = [tx_500, tx_500_lim[0], tx_500_lim[1]]
    zed[i] = z

    e_z =  sqrt(omega_m*(1.0 + z)^3 + omega_lam)
    yaxis = alog10(e_z^(4./3.) * ent_500)
    xaxis = alog10(tx_500/5.0)
    xup = alog10(tx_500_lim[1]/5.0) & xdown = alog10(tx_500_lim[0]/5.0)
    yup = alog10(e_z^(4./3.) * ent_500_lim[1]) 
    ydown = alog10(e_z^(4./3.) * ent_500_lim[0]) 
    xerr = 0.5*abs(xup - xdown) & yerr = 0.5*abs(yup - ydown)

    printf, 1, xaxis, xerr, yaxis, yerr, ey1ey2
endfor
close, 1
save, clu_names, ent_val, temp_val, zed, filename='ent_500.sav_part'
end

;Run BCES to get best fit
;part results
;2500
;B = 0.84(0.03) , A = 2.879(0.006)  Y|X
;B = 0.86(0.03) , A = 2.879(0.006)  Orth
;500
;B = 0.90(0.12) , A = 3.13(0.02)  Y|X
;B = 1.06(0.17) , A = 3.12(0.02)  Orth


pro make_figs
omega_m = 0.3 & omega_lam = 0.7

;2500 figure

restore, 'ent_2500.sav_part'
n = n_elements(zed)

e_z =  sqrt(omega_m*(1.0 + zed)^3 + omega_lam)
yaxis = e_z^(4./3.) * ent_val[0,*]
xaxis = temp_val[0,*]
xup = temp_val[2,*] & xdown = temp_val[1,*]
yup = e_z^(4./3.) * ent_val[2,*] 
ydown = e_z^(4./3.) * ent_val[1,*] 

set_plot, 'ps'
device, filename='ent_temp_2500.ps', /inches, xsize=5.3, /color
!p.multi=[0,1,1]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

plot, xaxis, yaxis, psym=4, /xlog, /ylog, xrange=[1,20], xstyle=1, $
  yrange=[100,3000], ystyle=1, xtitle='T!DX,2500!N (keV)', $
  ytitle='E(z)!U4/3!N K!D2500!N (keV cm!U2!N)' 
for i=0,n-1 do oplot, [xaxis[i],xaxis[i]], [ydown[i], yup[i]]
for i=0,n-1 do oplot, [xup[i],xdown[i]], [yaxis[i],yaxis[i]]

xfit = findgen(21)+1
yfit = 10.^(2.879) * (xfit/5)^(0.84)
oplot, xfit, yfit, thick=5.0, color=200

yfit = 10.^(2.879) * (xfit/5)^(0.86)
oplot, xfit, yfit, thick=5.0, color=50

device, /close

;500 figure

restore, 'ent_500.sav_part'
n = n_elements(zed)

e_z =  sqrt(omega_m*(1.0 + zed)^3 + omega_lam)
yaxis = e_z^(4./3.) * ent_val[0,*]
xaxis = temp_val[0,*]
xup = temp_val[2,*] & xdown = temp_val[1,*]
yup = e_z^(4./3.) * ent_val[2,*] 
ydown = e_z^(4./3.) * ent_val[1,*] 

device, filename='ent_temp_500.ps', /inches, xsize=5.3, /color
!p.multi=[0,1,1]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

plot, xaxis, yaxis, psym=4, /xlog, /ylog, xrange=[1,20], xstyle=1, $
  yrange=[300,5000], ystyle=1, xtitle='T!DX,500!N (keV)', $
  ytitle='E(z)!U4/3!N K!D500!N (keV cm!U2!N)' 
for i=0,n-1 do oplot, [xaxis[i],xaxis[i]], [ydown[i], yup[i]]
for i=0,n-1 do oplot, [xup[i],xdown[i]], [yaxis[i],yaxis[i]]
xyouts, 0.66, 950*5, '5000'
xyouts, 0.73, 950/2, '500'

xfit = findgen(21)+1
yfit = 10.^(3.13) * (xfit/5)^(0.90)
oplot, xfit, yfit, thick=5.0, color=200

yfit = 10.^(3.12) * (xfit/5)^(1.06)
oplot, xfit, yfit, thick=5.0, color=50

device, /close
set_plot, 'x'

end


pro calc_scat_k
omega_m = 0.3 & omega_lam = 0.7

;2500 scat

restore, 'ent_2500.sav_part'
n = n_elements(zed)

e_z =  sqrt(omega_m*(1.0 + zed)^3 + omega_lam)
yaxis = e_z^(4./3.) * ent_val[0,*]
xaxis = temp_val[0,*]
xup = temp_val[2,*] & xdown = temp_val[1,*]
yup = e_z^(4./3.) * ent_val[2,*] 
ydown = e_z^(4./3.) * ent_val[1,*] 
xerr = 0.5*abs(xup - xdown) & yerr = 0.5*abs(yup - ydown)
    
b = [2.879, 2.879]
m = [0.84, 0.86]

for i=0,1 do begin
    val = ((yaxis - (10.^b[i]*(xaxis/5.0)^m[i]))^2 - $
           (yerr)^2)/(yaxis)^2
    scat = sqrt(total(val)/(n-2.0))

    print, '2500 ', scat
endfor

;2500 scat

restore, 'ent_500.sav_part'
n = n_elements(zed)

e_z =  sqrt(omega_m*(1.0 + zed)^3 + omega_lam)
yaxis = e_z^(4./3.) * ent_val[0,*]
xaxis = temp_val[0,*]
xup = temp_val[2,*] & xdown = temp_val[1,*]
yup = e_z^(4./3.) * ent_val[2,*] 
ydown = e_z^(4./3.) * ent_val[1,*] 
xerr = 0.5*abs(xup - xdown) & yerr = 0.5*abs(yup - ydown)
    
b = [3.13, 3.12]
m = [0.90, 1.06]

for i=0,1 do begin
    val = ((yaxis - (10.^(b[i])*(xaxis/5.0)^m[i]))^2 - $
           (yerr)^2)/(yaxis)^2
    scat = sqrt(total(val)/(n-2.0))

    print, '500 ', scat
endfor
end
;;10% at r2500, ~30% at 500
