pro get_basic

files_all = ' '
readcol, 'ne_save/ne_files', files_all, format='A'
nc = n_elements(files_all)

openw, 1, 'clu_names'

for q = 0,nc-1 do begin
    restore, 'ne_save/'+files_all[q]
    ne_best_fit = parinfo.value
    x = where((status_all gt 0) and (status_all lt 5))
    if (x(0) ne -1) then ne_sample = results[x,*] else $
      print, name, 'No ne err fits'

    save, area, da, da_cm, dl, name, norm, normerr, $
      norm2, norm2err, n_m, ne_best_fit, ne_sample, rin, rin_cm, rin_mpc, $
      rout, rout_cm, rout_mpc, r_mpc, z, $
      filename=strcompress(name,/remove_all)+'.sav_t1'

    restore, 'temp_save/'+files_all[q]
    temp_best_fit = parinfo.value
    x = where((status_all gt 0) and (status_all lt 5))
    if (x(0) ne -1) then temp_sample = results[x,*] else $
      print, name, 'No temp err fits'
   
    save, abund, abunderr, calfile, temp_best_fit, temp_sample, temp, $
      temperr, filename=strcompress(name,/remove_all)+'.sav_t2'

    printf, 1, name
print, q
endfor
print, q

close, 1
end


pro get_basic_double

files_all = ['A1835+A1835aerr.sav', 'Perseus+A426err.sav',$
             'Sersic159-3+AS1101err.sav'] 
nc = n_elements(files_all)

for q = 0,nc-1 do begin
    restore, 'ne_save/'+files_all[q]
    ne_best_fit = parinfo.value
    x = where((status_all gt 0) and (status_all lt 5))
    if (x(0) ne -1) then ne_sample = results[x,*] else $
      print, name, 'No ne err fits'

    save, da, da_cm, dl, name, ne_best_fit, ne_sample, z, $
      filename=strcompress(name,/remove_all)+'.sav_t1'

    restore, 'temp_save/'+files_all[q]
    temp_best_fit = parinfo.value
    x = where((status_all gt 0) and (status_all lt 5))
    if (x(0) ne -1) then temp_sample = results[x,*] else $
      print, name, 'No temp err fits'

    ;calc area, norm2, norm2err, 
    area = !pi*((rout/60.)^2 - (rin/60.)^2)
    norm2 = norm*area
    norm2err = normerr*area 
   
    save, area, abund, abunderr, calfile, norm, normerr, norm2, norm2err, $
      n_m, rin, rin_cm, rin_mpc, rout, rout_cm, rout_mpc, r_mpc, $
      temp_best_fit, temp_sample, temp, temperr, $
      filename=strcompress(name,/remove_all)+'.sav_t2'

endfor

end

pro combine_basic
clu_names = ''
readcol, 'clu_names', clu_names, format='A'
n = n_elements(clu_names)

for q = 0,n-1 do begin
    restore, strcompress(clu_names[q],/remove_all)+'.sav_t2' 
    restore, strcompress(clu_names[q],/remove_all)+'.sav_t1'

    save, abund, abunderr, area, da, da_cm, dl, name, norm, normerr, $
      norm2, norm2err, n_m, ne_best_fit, ne_sample, rin, rin_cm, rin_mpc, $
      rout, rout_cm, rout_mpc, r_mpc, z, calfile, temp_best_fit, $
      temp_sample, temp, temperr, $
      filename=strcompress(name,/remove_all)+'_basic.sav'
endfor
end


;;;Do mass, gas mass, temperature, entropy as fxn of radius &
;;;calculate the errors for each
pro radial_dists
clu_names = ''
readcol, 'clu_names', clu_names, format='A'
n = n_elements(clu_names)

for q = 0,n-1 do begin
    restore, strcompress(clu_names[q],/remove_all)+'_basic.sav'

    n_temps = n_elements(temp_sample[*,0])
    n_nes = n_elements(ne_sample[*,0])

    tempdist = fltarr(n_m,3)
    tempdist[*,0] = temp_dist(r_mpc, temp_best_fit)

    td_t = fltarr(n_m,n_temps)
    for i=0,n_temps-1 do td_t[*,i] = temp_dist(r_mpc, temp_sample[i,*])
    td_r = fltarr(n_m,2)
    for i=0,n_m-1 do td_r[i,*] = get_limit(tempdist[i,0], td_t[i,*], 0.8413)
    tempdist[*,1] = td_r[*,0]
    tempdist[*,2] = td_r[*,1]

    nedist = fltarr(n_m,3)
    nedist[*,0] = ne_dist(r_mpc, ne_best_fit)

    ned_t = fltarr(n_m,n_nes)
    for i=0,n_nes-1 do ned_t[*,i] = ne_dist(r_mpc, ne_sample[i,*])
    ned_r = fltarr(n_m,2)
    for i=0,n_m-1 do ned_r[i,*] = get_limit(nedist[i,0], ned_t[i,*], 0.8413)
    nedist[*,1] = ned_r[*,0]
    nedist[*,2] = ned_r[*,1]

    mgdist = fltarr(n_m,3)
    mgdist[*,0] = gas_mass(r_mpc, ne_best_fit)

    mgd_t = fltarr(n_m,n_nes)
    for i=0,n_nes-1 do mgd_t[*,i] = gas_mass(r_mpc, ne_sample[i,*])
    mgd_r = fltarr(n_m,2)
    for i=0,n_m-1 do mgd_r[i,*] = get_limit(mgdist[i,0], mgd_t[i,*], 0.8413)
    mgdist[*,1] = mgd_r[*,0]
    mgdist[*,2] = mgd_r[*,1]

    massdist = fltarr(n_m,3)
    massdist[*,0] = total_mass(r_mpc, ne_best_fit, temp_best_fit)

    md_t = fltarr(n_m,n_nes,n_temps)
    for i=0,n_nes-1 do $
      for j=0,n_temps-1 do $
        md_t[*,i] = total_mass(r_mpc, ne_sample[i,*], temp_sample[j,*])
    md_r = fltarr(n_m,2)
    for i=0,n_m-1 do md_r[i,*] = get_limit(massdist[i,0], md_t[i,*,*], 0.8413)
    massdist[*,1] = md_r[*,0]
    massdist[*,2] = md_r[*,1]

    entdist = fltarr(n_m,3)
    entdist[*,0] = entropy_dist(r_mpc, ne_best_fit, temp_best_fit)

    entd_t = fltarr(n_m,n_nes,n_temps)
    for i=0,n_nes-1 do $
      for j=0,n_temps-1 do $
        entd_t[*,i] = entropy_dist(r_mpc, ne_sample[i,*], temp_sample[j,*])
    entd_r = fltarr(n_m,2)
    for i=0,n_m-1 do entd_r[i,*] = $
      get_limit(entdist[i,0], entd_t[i,*,*], 0.8413)
    entdist[*,1] = entd_r[*,0]
    entdist[*,2] = entd_r[*,1]

    save, r_mpc, tempdist, nedist, mgdist, massdist, entdist, name, n_m, $
      filename=strcompress(name,/remove_all)+'_dists.sav'

endfor
end


;Next step - calculate r2500 & r500 for each.  Output cluster names to
;            files: clu_2500 & clu_500

pro calc_r_vals
clu_names = ''
readcol, 'clu_names', clu_names, format='A'
n = n_elements(clu_names)

openw, 1, 'clu_2500'
openw, 2, 'clu_500'

set_plot, 'ps'
device, filename='radii_calc.ps', /inches, ysize=9.0, yoffset=1.0
!p.multi=[0,1,2]

for q=0,n-1 do begin
    r2500 = -999 & r500 = -999
    restore, strcompress(clu_names[q],/remove_all)+'_basic.sav'
    r = [r_mpc,max(rout_mpc)]

;delta = 2500
    del = 2500.
    res = mass_comp(r, ne_best_fit, temp_best_fit, del, z)
    plot, r, res, psym=4, title=name
    oplot, [0,100], [0,0]

if (strcompress(name,/remove_all) eq 'A2163') then $
  print, 'skipping A2163 r2500' $
  else begin
    x = where(res lt 0.0)
    if (x[0] ne -1) then begin
        if (n_elements(x) gt 1) then begin
            x2 = where(res(x) gt res(x+1))
            x = x[x2]
        endif
        r_test = findgen(51)/50.*(r[x(0)]-r[x(0)-1]) + r[x(0)-1]
        res2 = mass_comp(r_test, ne_best_fit, temp_best_fit, del, z)
        
        y = where(res2 lt 0.0)
        r2500 = r_test[y(0)] - res2[y(0)]*$
          (r_test[y(0)] - r_test[y(0)-1])/(res2[y(0)] - res2[y(0)-1])

        printf, 1, name
    endif
endelse

;delta = 500
    del = 500.
    res = mass_comp(r, ne_best_fit, temp_best_fit, del, z)
    plot, r, res, psym=4, title=name
    oplot, [0,100], [0,0]

    x = where(res lt 0.0)
    if (x[0] ne -1) then begin
        r_test = findgen(51)/50.*(r[x(0)]-r[x(0)-1]) + r[x(0)-1]
        res2 = mass_comp(r_test, ne_best_fit, temp_best_fit, del, z)

        y = where(res2 lt 0.0)
        r500 = r_test[y(0)] - res2[y(0)]*$
          (r_test[y(0)] - r_test[y(0)-1])/(res2[y(0)] - res2[y(0)-1])

        printf, 2, name
    endif
    
    save, r500, r2500, name, $
      filename=strcompress(name,/remove_all)+'_derive.sav'

endfor

close, 1
close, 2
device, /close
set_plot, 'x'

end

;Check r500 & r2500 values
pro get_rvals

clu_names = ''
readcol, 'clu_names', clu_names, format='A'
n = n_elements(clu_names)

r_2500 = fltarr(n)
r_500 = fltarr(n)
names = strarr(n)
for q=0,n-1 do begin
    restore, strcompress(clu_names[q],/remove_all)+'_derive.sav'

    r_2500[q] = r2500
    r_500[q] = r500
    names[q] = name
endfor

stop
end
;mean r2500/r500 = 0.45 --> rin = 0.33r2500


;For r500 systems -- calc r_err, mass (total & gas), mass errs, temp
;                    (0.15-1 r500), temp errs, yx, yx errs, fg, fg
;                    errs

pro calc_500
common int_issues, rmax, ne_pars, temp_pars, abund, rin_mpc, rout_mpc

;r500
files = 'a'
readcol, 'clu_500', files, format='A'
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

        r500_err = fltarr(n_nes, n_temps)
        r = [r_mpc,max(rout_mpc)]

 ;delta = 500
        del = 500.
        for j=0,n_nes-1 do begin
            for k=0,n_temps-1 do begin
                res = mass_comp(r, ne_sample[j,*], temp_sample[k,*], del, z)
                x = where(res lt 0.0)
                if (x[0] eq -1) then r500_err[j,k] = max(r) $
                else if ((n_elements(x) eq 1) and (x[0] le 2)) then $
                  r500_err[j,k] = max(r) $
                else begin
                    if (x[0] le 2) then x = x[1:*]
                    r_test = findgen(51)/50.*(r[x(0)]-r[x(0)-1]) + r[x(0)-1]
                    res2 = mass_comp(r_test, ne_sample[j,*], $
                                 temp_sample[k,*], del, z)
                    y = where(res2 lt 0.0)
                    rtemp = r_test[y(0)] - res2[y(0)]*(r_test[y(0)] $
                        - r_test[y(0)-1])/(res2[y(0)] - res2[y(0)-1])
                
                    r500_err[j,k] = rtemp
                endelse
            endfor
        endfor
        
        r500_lim = get_limit(r500, r500_err[*,*], 0.8413)
        res = [r500, r500_lim[0], r500_lim[1]]
    
        mass_tot_500 = total_mass(r500, ne_best_fit, temp_best_fit)
        mass_gas_500 = gas_mass(r500, ne_best_fit)

        rmax = max(rout_mpc)
        ne_pars = ne_best_fit
        temp_pars = temp_best_fit
        tx_500 = calc_temp_spec_v2('xmm/mos', calfile, $
                                   0.15*r500, r500)

        mgasb = gas_mass(res, ne_best_fit)
        mtotb = total_mass(res, ne_best_fit, temp_best_fit)
        txb = fltarr(3)
        for j=0,2 do txb[j] = calc_temp_spec_v2('xmm/mos', calfile, $
                                   0.15*res[j], res[j])

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
                                   0.15*res[j], res[j])
                endfor
            endfor
        endfor

        yx_500 = tx_500 * 10.^(mass_gas_500)
        fg_500 = 10.^(mass_gas_500 - mass_tot_500)

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

        mass_tot_500_lim = fltarr(2)
        mass_gas_500_lim = fltarr(2)
        tx_500_lim = fltarr(2)
        yx_500_lim = fltarr(2)
        fg_500_lim = fltarr(2)

        x = where(mtot_errtemp[*,0] gt 0)
        mass_tot_500_lim[0] = min(mtot_errtemp[x,0])
        mass_tot_500_lim[1] = max(mtot_errtemp[*,1])
   
        x = where(mgas_errtemp[*,0] gt 0)
        mass_gas_500_lim[0] = min(mgas_errtemp[x,0])
        mass_gas_500_lim[1] = max(mgas_errtemp[*,1])

        x = where(tx_errtemp[*,0] gt 0)
        tx_500_lim[0] = min(tx_errtemp[x,0])
        tx_500_lim[1] = max(tx_errtemp[*,1])

        x = where(yx_errtemp[*,0] gt 0)
        yx_500_lim[0] = min(yx_errtemp[x,0])
        yx_500_lim[1] = max(yx_errtemp[*,1])

        x = where(fg_errtemp[*,0] gt 0)
        fg_500_lim[0] = min(fg_errtemp[x,0])
        fg_500_lim[1] = max(fg_errtemp[*,1])


        save, r500, name, r500_err, r500_lim, mass_tot_500, mass_gas_500, $
          tx_500, yx_500, fg_500, mass_tot_err, mass_gas_err, tx_err, $
          yx_err, fg_err, mass_tot_500_lim, mass_gas_500_lim, tx_500_lim, $
          yx_500_lim, fg_500_lim, $
          filename=strcompress(name,/remove_all)+'_500derive.sav'
    endelse

endfor

end
