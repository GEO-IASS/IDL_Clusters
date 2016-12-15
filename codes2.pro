pro calc_2500
common int_issues, rmax, ne_pars, temp_pars, abund, rin_mpc, rout_mpc

;r2500
files = 'a'
readcol, 'clu_2500', files, format='A'
n_clu = n_elements(files)

for i=0,n_clu-1 do begin
    restore, strcompress(files[i],/remove_all)+'_basic.sav'
    restore, strcompress(files[i],/remove_all)+'_derive.sav'

    if ((temp_sample[0,0] eq -999) or (ne_sample[0,0] eq -999)) then begin
        print, 'Error analysis not possible, missing error information.' 
        radii_err = -999
    endif else begin
        n_temps = n_elements(temp_sample[*,0])
        n_nes = n_elements(ne_sample[*,0])

        r2500_err = fltarr(n_nes, n_temps)
        r = [r_mpc,max(rout_mpc)]

 ;delta = 2500
        del = 2500.
        for j=0,n_nes-1 do begin
            for k=0,n_temps-1 do begin
                res = mass_comp(r, ne_sample[j,*], temp_sample[k,*], del, z)
                x = where(res lt 0.0)
                if (x[0] eq -1) then r2500_err[j,k] = max(r) $
                else if ((n_elements(x) eq 1) and (x[0] le 2)) then $
                  r2500_err[j,k] = max(r) $
                else begin
                    if (x[0] le 2) then x = x[1:*]
                    r_test = findgen(51)/50.*(r[x(0)]-r[x(0)-1]) + r[x(0)-1]
                    res2 = mass_comp(r_test, ne_sample[j,*], $
                                 temp_sample[k,*], del, z)
                    y = where(res2 lt 0.0)
                    rtemp = r_test[y(0)] - res2[y(0)]*(r_test[y(0)] $
                        - r_test[y(0)-1])/(res2[y(0)] - res2[y(0)-1])
                
                    r2500_err[j,k] = rtemp
                endelse
            endfor
        endfor
        
        r2500_lim = get_limit(r2500, r2500_err[*,*], 0.8413)
        res = [r2500, r2500_lim[0], r2500_lim[1]]
    
        mass_tot_2500 = total_mass(r2500, ne_best_fit, temp_best_fit)
        mass_gas_2500 = gas_mass(r2500, ne_best_fit)

        rmax = max(rout_mpc)
        ne_pars = ne_best_fit
        temp_pars = temp_best_fit
        tx_2500 = calc_temp_spec_v2('xmm/mos', calfile, $
                                   0.33*r2500, r2500)

        mgasb = gas_mass(res, ne_best_fit)
        mtotb = total_mass(res, ne_best_fit, temp_best_fit)
        txb = fltarr(3)
        for j=0,2 do txb[j] = calc_temp_spec_v2('xmm/mos', calfile, $
                                   0.33*res[j], res[j])

        mass_tot_err = fltarr(3,n_nes,n_temps)
        mass_gas_err = fltarr(3,n_nes,n_temps)
        tx_err = fltarr(3, n_nes, n_temps)
        for j=0,2 do begin
            for k=0,n_nes-1 do begin
                ne_pars = ne_sample[k,*]
                for l=0,n_temps-1 do begin
                    temp_pars = temp_sample[l,*]

                    mass_tot_err[j,k,l] = total_mass(res[j], ne_sample[k,*], $
                                                     temp_sample[l,*])
                    mass_gas_err[j,k,l] = gas_mass(res[j], ne_sample[k,*])

                    tx_err[j,k,l] = calc_temp_spec_v2('xmm/mos', calfile, $
                                   0.33*res[j], res[j])
                endfor
            endfor
        endfor

        yx_2500 = tx_2500 * 10.^(mass_gas_2500)
        fg_2500 = 10.^(mass_gas_2500 - mass_tot_2500)

        yxb = txb * 10.^(mgasb)
        fgb = 10.^(mgasb - mtotb)

        yx_err = tx_err * 10.^(mass_gas_err)
        fg_err = 10.^(mass_gas_err - mass_tot_err)
        
;Now get limits
        mtot_errtemp = fltarr(3,2)
        mgas_errtemp = fltarr(3,2)
        tx_errtemp = fltarr(3,2)
        yx_errtemp = fltarr(3,2)
        fg_errtemp = fltarr(3,2)
        for j=0,2 do begin
            mtot_errtemp[j,*] = get_limit(mtotb[j], mass_tot_err[j,*,*], $
                                          0.8413)
            mgas_errtemp[j,*] = get_limit(mgasb[j], mass_gas_err[j,*,*], $
                                          0.8413)
            tx_errtemp[j,*] = get_limit(txb[j], tx_err[j,*,*], 0.8413)
            yx_errtemp[j,*] = get_limit(yxb[j], yx_err[j,*,*], 0.8413)
            fg_errtemp[j,*] = get_limit(fgb[j], fg_err[j,*,*], 0.8413)
        endfor

        mass_tot_2500_lim = fltarr(2)
        mass_gas_2500_lim = fltarr(2)
        tx_2500_lim = fltarr(2)
        yx_2500_lim = fltarr(2)
        fg_2500_lim = fltarr(2)

        x = where(mtot_errtemp[*,0] gt 0)
        mass_tot_2500_lim[0] = min(mtot_errtemp[x,0])
        mass_tot_2500_lim[1] = max(mtot_errtemp[*,1])
   
        x = where(mgas_errtemp[*,0] gt 0)
        mass_gas_2500_lim[0] = min(mgas_errtemp[x,0])
        mass_gas_2500_lim[1] = max(mgas_errtemp[*,1])

        x = where(tx_errtemp[*,0] gt 0)
        tx_2500_lim[0] = min(tx_errtemp[x,0])
        tx_2500_lim[1] = max(tx_errtemp[*,1])

        x = where(yx_errtemp[*,0] gt 0)
        yx_2500_lim[0] = min(yx_errtemp[x,0])
        yx_2500_lim[1] = max(yx_errtemp[*,1])

        x = where(fg_errtemp[*,0] gt 0)
        fg_2500_lim[0] = min(fg_errtemp[x,0])
        fg_2500_lim[1] = max(fg_errtemp[*,1])


        save, r2500, name, r2500_err, r2500_lim, mass_tot_2500, $
          mass_gas_2500, $
          tx_2500, yx_2500, fg_2500, mass_tot_err, mass_gas_err, tx_err, $
          yx_err, fg_err, mass_tot_2500_lim, mass_gas_2500_lim, tx_2500_lim, $
          yx_2500_lim, fg_2500_lim, $
          filename=strcompress(name,/remove_all)+'_2500derive.sav'
    endelse

endfor

end


;Add M500, Mgas500, tx500, yx500, fg500 & same for 2500 in derive

pro calc_derived
common int_issues, rmax, ne_pars, temp_pars, abund, rin_mpc, rout_mpc

files = 'a'
readcol, 'clu_names', files, format='A'
n_clu = n_elements(files)

for i=0,n_clu-1 do begin
    mass_tot_2500 = -999 & mass_tot_500 = -999
    mass_gas_2500 = -999 & mass_gas_500 = -999
    tx_2500 = -999 & tx_500 = -999
    yx_2500 = -999 & yx_500 = -999
    fg_2500 = -999 & fg_500 = -999

    restore, strcompress(files[i],/remove_all)+'_basic.sav'
    restore, strcompress(files[i],/remove_all)+'_derive.sav'

    if (r2500 gt 0) then begin
        mass_tot_2500 = total_mass(r2500, ne_best_fit, temp_best_fit)
        mass_gas_2500 = gas_mass(r2500, ne_best_fit)

        rmax = max(rout_mpc)
        ne_pars = ne_best_fit
        temp_pars = temp_best_fit
        tx_2500 = calc_temp_spec_v2('xmm/mos', calfile, 0.33*r2500, r2500)
        yx_2500 = tx_2500 * 10.^(mass_gas_2500)
        fg_2500 = 10.^(mass_gas_2500 - mass_tot_2500)
    endif

    if (r500 gt 0) then begin
        mass_tot_500 = total_mass(r500, ne_best_fit, temp_best_fit)
        mass_gas_500 = gas_mass(r500, ne_best_fit)

        rmax = max(rout_mpc)
        ne_pars = ne_best_fit
        temp_pars = temp_best_fit
        tx_500 = calc_temp_spec_v2('xmm/mos', calfile, 0.15*r500, r500)
        yx_500 = tx_500 * 10.^(mass_gas_500)
        fg_500 = 10.^(mass_gas_500 - mass_tot_500)
    endif

    save, r500, r2500, name, mass_tot_2500, mass_gas_2500, $
      tx_2500, yx_2500, fg_2500, mass_tot_500, mass_gas_500, $
      tx_500, yx_500, fg_500, $
      filename=strcompress(name,/remove_all)+'_derive.sav'
endfor

end


;;Do for scaled projected temperature profiles.
;;need temperature result
;;tx, r500

pro make_tempradprof

files = 'a'
readcol, 'clu_500', files, format='A'
n_clu = n_elements(files)

set_plot,'ps'
device, filename='radprof_fig.ps'
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0 
!p.charthick=2.0

restore, strcompress(files[0],/remove_all)+'_basic.sav'
;get temp, temperr, r_mpc, rin_mpc, rout_mpc
restore, strcompress(files[0],/remove_all)+'_derive.sav'
;get tx, r500

plot, r_mpc/r500, temp[*]/tx_500, psym=3, xrange=[0,2.0], yrange=[0,2.0], $
  ytitle='T / T!BX!N', xtitle='Radius (r!B500!N)', xmargin=[8,3]

for j=0,n_m-1 do oplot, [r_mpc[j],r_mpc[j]]/r500, $
  (temp[j] + [temperr[j],-temperr[j]])/tx_500
for j=0,n_m-1 do oplot, [rin_mpc[j],rout_mpc[j]]/r500, $
  [temp[j],temp[j]]/tx_500

for i=1,n_clu-1 do begin
    restore, strcompress(files[i],/remove_all)+'_basic.sav'
;get temp, temperr, r_mpc, rin_mpc, rout_mpc
    restore, strcompress(files[i],/remove_all)+'_derive.sav'
;get tx, r500

    oplot, r_mpc/r500, temp[*]/tx_500, psym=3
    for j=0,n_m-1 do oplot, [r_mpc[j],r_mpc[j]]/r500, $
      (temp[j] + [temperr[j],-temperr[j]])/tx_500
    for j=0,n_m-1 do oplot, [rin_mpc[j],rout_mpc[j]]/r500, $
      [temp[j],temp[j]]/tx_500

endfor
device, /close
set_plot, 'x'
end


pro make_avgtdist_fig
file_name = 'a'
readcol, 'clu_500', file_name, format='A'
n_file = n_elements(file_name)

rbin_low = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3]
rbin_hi = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,2.0]
nbins = n_elements(rbin_low)

t_sum = fltarr(nbins)
t_wt = fltarr(nbins)
bin_tvals = fltarr(nbins,500)
bin_tvals_cnt = fltarr(nbins)

rbin_mid = 0.5*(rbin_low+rbin_hi)

for k=1,n_file-1 do begin
    restore, strcompress(file_name[k],/remove_all)+'_basic.sav'
;get temp, temperr, r_mpc, rin_mpc, rout_mpc
    restore, strcompress(file_name[k],/remove_all)+'_derive.sav'
;get tx, r500

    for i=0,nbins-1 do begin
    ;type 1 -- fully contained
        x = where((rin_mpc/r500 ge rbin_low[i]) and $
                  (rout_mpc/r500 le rbin_hi[i]),c)
    
        if (c gt 0) then begin
            t_sum[i] = t_sum[i] + $
              total(temp[x]/tx_500/temperr[x]/temperr[x])
            t_wt[i] = t_wt[i] + total(1.0/temperr[x]/temperr[x])
            bin_tvals_cnt[i] = bin_tvals_cnt[i] + c
        endif

    ;type 2 -- in<binlow, out>binlow<binhi
        x = where((rin_mpc/r500 lt rbin_low[i]) and $
                  ((rout_mpc/r500 le rbin_hi[i]) and $
                   (rout_mpc/r500 gt rbin_low[i])),c)
        if (c gt 0) then begin
            for j=0,c-1 do begin
                wt = (rout_mpc[x[j]]/r500 - rbin_low[i])/$
                  ((rout_mpc[x[j]]-rin_mpc[x[j]])/r500)
                
                t_sum[i] = t_sum[i] + $
                  temp[x[j]]*wt/tx_500/temperr[x[j]]/temperr[x[j]]
                t_wt[i] = t_wt[i] + wt/temperr[x[j]]/temperr[x[j]]
                bin_tvals_cnt[i] = bin_tvals_cnt[i] + wt
            endfor
        endif

    ;type 3 -- out>binout, in>binlow<binhi
        x = where((rout_mpc/r500 gt rbin_hi[i]) and $
                  ((rin_mpc/r500 lt rbin_hi[i]) and $
                   (rin_mpc/r500 ge rbin_low[i])),c)
        if (c gt 0) then begin
            for j=0,c-1 do begin
                wt = (rbin_hi[i]-rin_mpc[x[j]]/r500)/$
                  ((rout_mpc[x[j]]-rin_mpc[x[j]])/r500)

                t_sum[i] = t_sum[i] + $
                  temp[x[j]]*wt/tx_500/temperr[x[j]]/temperr[x[j]]
                t_wt[i] = t_wt[i] + wt/temperr[x[j]]/temperr[x[j]]
                bin_tvals_cnt[i] = bin_tvals_cnt[i] + wt
            endfor
        endif

    ;type 4 -- in<binlow, out>binhi
        x = where((rout_mpc/r500 gt rbin_hi[i]) and $
                  (rin_mpc/r500 lt rbin_low[i]),c)
        if (c gt 0) then begin
            for j=0,c-1 do begin
                wt = (rbin_hi[i]-rbin_low[i])/$
                  ((rout_mpc[x[j]]-rin_mpc[x[j]])/r500)

                t_sum[i] = t_sum[i] + $
                  temp[x[j]]*wt/tx_500/temperr[x[j]]/temperr[x[j]]
                t_wt[i] = t_wt[i] + wt/temperr[x[j]]/temperr[x[j]]
                bin_tvals_cnt[i] = bin_tvals_cnt[i] + wt
            endfor
        endif

    ;type 5 -- in&out<binlow - do nothing
    ;type 6 -- in&out>binhi - do nothing

    endfor
endfor

t_avg = t_sum/t_wt
t_avg_err = sqrt(1.0/t_wt)

;;determine the spread in values.
vari = fltarr(nbins)
for i=0,nbins-1 do begin
    vals = 0
    for k=1,n_file-1 do begin
        restore, strcompress(file_name[k],/remove_all)+'_basic.sav'
;get temp, temperr, r_mpc, rin_mpc, rout_mpc
        restore, strcompress(file_name[k],/remove_all)+'_derive.sav'
;get tx, r500

    ;type 1 -- fully contained
        x = where((rin_mpc/r500 ge rbin_low[i]) and $
                  (rout_mpc/r500 le rbin_hi[i]),c)
        if (c gt 0) then vals = [vals, temp[x]/tx_500]

    ;type 2 -- in<binlow, out>binlow<binhi
        x = where((rin_mpc/r500 lt rbin_low[i]) and $
                  ((rout_mpc/r500 le rbin_hi[i]) and $
                   (rout_mpc/r500 gt rbin_low[i])),c)
        if (c gt 0) then vals = [vals, temp[x]/tx_500]
       
    ;type 3 -- out>binout, in>binlow<binhi
        x = where((rout_mpc/r500 gt rbin_hi[i]) and $
                  ((rin_mpc/r500 lt rbin_hi[i]) and $
                   (rin_mpc/r500 ge rbin_low[i])),c)
        if (c gt 0) then vals = [vals, temp[x]/tx_500]

    ;type 4 -- in<binlow, out>binhi
        x = where((rout_mpc/r500 gt rbin_hi[i]) and $
                  (rin_mpc/r500 lt rbin_low[i]),c)
        if (c gt 0) then vals = [vals, temp[x]/tx_500]

    ;type 5 -- in&out<binlow - do nothing
    ;type 6 -- in&out>binhi - do nothing

    endfor
    
    vari[i] = sqrt(total((vals[1:*] - t_avg[i])^2)/(n_elements(vals)-2))
endfor


set_plot, 'ps'
device, filename='avgtdist_fig.ps'
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0 
!p.charthick=2.0

plot, rbin_mid, t_avg, psym=5, yrange=[0,1.5], xtitle='r/r_500', ytitle='T/T_X'
for j=0,nbins-1 do oplot, [rbin_mid[j],rbin_mid[j]], $
  t_avg[j]+[-t_avg_err[j],t_avg_err[j]]
for j=0,nbins-1 do oplot, [rbin_low[j],rbin_hi[j]], [t_avg[j],t_avg[j]]
oplot, rbin_mid, t_avg+vari, linestyle=1.0
oplot, rbin_mid, t_avg-vari, linestyle=1.0

;Fit linear model r>0.3
x = where(rbin_mid gt 0.3)

Result = LINFIT(rbin_mid[x], t_avg[x], MEASURE_ERRORS=t_avg_err[x], SIGMA=sigma, YFIT=yfit)

print, result, sigma
 
oplot, rbin_mid[x], yfit

device, /close
set_plot, 'x'

;t/tx = 1.35(0.03) - 0.78(0.05)(r/r500)


end

        
