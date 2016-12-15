pro test_fig

files ='a'
readcol, 'clu_names', files, format='A'
n_clu = n_elements(files)

set_plot, 'ps'
device, filename='abundvstemp.ps', /color
loadct,6

;For first cluster
restore, strcompress(files[0],/remove_all)+'_basic.sav'
x = where(abunderr eq 0.0,c)
if (x[0] ne -1) then for j=0,c-1 do abunderr(x[j]) = abunderr(x[j]-1)

plot, temp, abund, psym=4, xrange=[0,25], yrange=[0,2], $
  xtitle='Temperature (keV)', ytitle='Abundance'
for j=0,n_m-1 do oplot, [temp[j],temp[j]], abund[j]+[-abunderr[j],abunderr[j]]
for j=0,n_m-1 do oplot, temp[j]+[-temperr[j],temperr[j]], [abund[j],abund[j]]

for i=1,n_clu-1 do begin
    restore, strcompress(files[i],/remove_all)+'_basic.sav'
    x = where(abunderr eq 0.0,c)
    if (x[0] ne -1) then for j=0,c-1 do abunderr(x[j]) = abunderr(x[j]-1)

    oplot, temp, abund, psym=4
    for j=0,n_m-1 do oplot, [temp[j],temp[j]], abund[j]+[-abunderr[j],abunderr[j]]
    for j=0,n_m-1 do oplot, temp[j]+[-temperr[j],temperr[j]], [abund[j],abund[j]]

endfor

device, /close
set_plot, 'x'

end


;;bin by temperature.  use mean weight.
pro binned_tempabund
files ='a'
readcol, 'clu_names', files, format='A'
n_clu = n_elements(files)

tbin_low = [0.5,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,12.5]
tbin_hi = [1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,12.5,25]
nbins = n_elements(tbin_low)

a_sum = fltarr(nbins)
a_wt = fltarr(nbins)
bin_avals = fltarr(nbins,500)
bin_avals_cnt = fltarr(nbins)

tbin_mid = 0.5*(tbin_low+tbin_hi)

for k=1,n_clu-1 do begin
    restore, strcompress(files[k],/remove_all)+'_basic.sav'
    x = where(abunderr eq 0.0,c)
    if (x[0] ne -1) then for j=0,c-1 do abunderr(x[j]) = abunderr(x[j]-1)

    tin = temp - temperr
    tout = temp + temperr

    for i=0,nbins-1 do begin
    ;type 1 -- fully contained
        x = where((tin ge tbin_low[i]) and (tout le tbin_hi[i]),c)
    
        if (c gt 0) then begin
            a_sum[i] = a_sum[i] + $
              total(abund[x]/abunderr[x]/abunderr[x])
            a_wt[i] = a_wt[i] + total(1.0/abunderr[x]/abunderr[x])
            bin_avals_cnt[i] = bin_avals_cnt[i] + c
        endif

    ;type 2 -- in<binlow, out>binlow<binhi
        x = where((tin lt tbin_low[i]) and ((tout le tbin_hi[i]) and $
                   (tout gt tbin_low[i])),c)
        if (c gt 0) then begin
            for j=0,c-1 do begin
                wt = (tout[x[j]] - tbin_low[i])/(tout[x[j]]-tin[x[j]])
                
                a_sum[i] = a_sum[i] + $
                  abund[x[j]]*wt/abunderr[x[j]]/abunderr[x[j]]
                a_wt[i] = a_wt[i] + wt/abunderr[x[j]]/abunderr[x[j]]
                bin_avals_cnt[i] = bin_avals_cnt[i] + wt
            endfor
        endif

    ;type 3 -- out>binout, in>binlow<binhi
        x = where((tout gt tbin_hi[i]) and ((tin lt tbin_hi[i]) and $
                   (tin ge tbin_low[i])),c)
        if (c gt 0) then begin
            for j=0,c-1 do begin
                wt = (tbin_hi[i]-tin[x[j]])/(tout[x[j]]-tin[x[j]])

                a_sum[i] = a_sum[i] + $
                  abund[x[j]]*wt/abunderr[x[j]]/abunderr[x[j]]
                a_wt[i] = a_wt[i] + wt/abunderr[x[j]]/abunderr[x[j]]
                bin_avals_cnt[i] = bin_avals_cnt[i] + wt
            endfor
        endif

    ;type 4 -- in<binlow, out>binhi
        x = where((tout gt tbin_hi[i]) and (tin lt tbin_low[i]),c)
        if (c gt 0) then begin
            for j=0,c-1 do begin
                wt = (tbin_hi[i]-tbin_low[i])/(tout[x[j]]-tin[x[j]])

                a_sum[i] = a_sum[i] + $
                  abund[x[j]]*wt/abunderr[x[j]]/abunderr[x[j]]
                a_wt[i] = a_wt[i] + wt/abunderr[x[j]]/abunderr[x[j]]
                bin_avals_cnt[i] = bin_avals_cnt[i] + wt
            endfor
        endif

    ;type 5 -- in&out<binlow - do nothing
    ;type 6 -- in&out>binhi - do nothing

    endfor
endfor
;print, bin_avals_cnt
a_avg = a_sum/a_wt
a_avg_err = sqrt(1.0/a_wt)

 ; determine the spread in values.
vari = fltarr(nbins)
for i=0,nbins-1 do begin
    vals = 0
    for k=1,n_clu-1 do begin
        restore, strcompress(files[k],/remove_all)+'_basic.sav'
        tin = temp - temperr
        tout = temp + temperr

     ;type 1 -- fully contained
         x = where((tin ge tbin_low[i]) and (tout le tbin_hi[i]),c)
         if (c gt 0) then vals = [vals, abund[x]]
;print, c
     ;type 2 -- in<binlow, out>binlow<binhi
         x = where((tin lt tbin_low[i]) and ((tout le tbin_hi[i]) and $
                    (tout gt tbin_low[i])),c)
         if (c gt 0) then vals = [vals, abund[x]]
;print, c       
     ;type 3 -- out>binout, in>binlow<binhi
         x = where((tout gt tbin_hi[i]) and ((tin lt tbin_hi[i]) and $
                    (tin ge tbin_low[i])),c)
         if (c gt 0) then vals = [vals, abund[x]]
;print, c
     ;type 4 -- in<binlow, out>binhi
         x = where((tout gt tbin_hi[i]) and (tin lt tbin_low[i]),c)
         if (c gt 0) then vals = [vals, abund[x]]
;print, c
     ;type 5 -- in&out<binlow - do nothing
     ;type 6 -- in&out>binhi - do nothing

     endfor
    
     vari[i] = sqrt(total((vals[1:*] - a_avg[i])^2)/(n_elements(vals)-2))
 endfor


set_plot, 'ps'
device, filename='abund_temp_bin.ps'
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0 
!p.charthick=2.0

sol = string(110B)
plot, tbin_mid, a_avg, psym=4, yrange=[0,1.0], xtitle='Temperature (keV)', $
  ytitle='Abundance (Z!D!9'+sol+'!X!N)'
for j=0,nbins-1 do oplot, [tbin_mid[j],tbin_mid[j]], $
  a_avg[j]+[-a_avg_err[j],a_avg_err[j]]
for j=0,nbins-1 do oplot, [tbin_low[j],tbin_hi[j]], [a_avg[j],a_avg[j]]
oplot, tbin_mid, a_avg+vari, linestyle=1.0
oplot, tbin_mid, a_avg-vari, linestyle=1.0

device, /close
set_plot, 'x'

end


function int_ner2, X
common int_fun, par
r = X
a = par

cusp = (r/a[1])^(-a[3])
brk = 1.0/(1.0 + r^3/a[4]^3)^(a[5]/3.) 
king = 1.0/(1.0 + r^2/a[1]^2)^(3*a[2] - a[3]/2.)
k2 = 1.0/(1.0 + r^2/a[7]^2)^(3*a[8])

ne2_fxn = a[0]*a[0]*cusp*brk*king + a[6]*a[6]*k2

ne_fxn = sqrt(ne2_fxn)

return, r*r*ne_fxn
end


function gas_m_int, rin, rout, a
common int_fun, par

par = a

nr = n_elements(rin)
int = fltarr(nr)
for i=0, nr-1 do begin
    if (rin[i] eq 0) then rin[i] = 1e-4
    int[i] = qromb('int_ner2', rin[i], rout[i], jmax=25)
endfor

return, int
end


function int_ne2r2, X
common int_fun, par
r = X
a = par

cusp = (r/a[1])^(-a[3])
brk = 1.0/(1.0 + r^3/a[4]^3)^(a[5]/3.) 
king = 1.0/(1.0 + r^2/a[1]^2)^(3*a[2] - a[3]/2.)
k2 = 1.0/(1.0 + r^2/a[7]^2)^(3*a[8])

ne2_fxn = a[0]*a[0]*cusp*brk*king + a[6]*a[6]*k2

return, r*r*ne2_fxn
end


function em_int, rin, rout, a
common int_fun, par

par = a

nr = n_elements(rin)
int = fltarr(nr)
for i=0, nr-1 do begin
    if (rin[i] eq 0) then rin[i] = 1e-4
    int[i] = qromb('int_ne2r2', rin[i], rout[i], jmax=25)
endfor

return, int
end


pro calc_weight_abund
;mass-weighted p
;em-weighted p^2
;int(wAdV)/int(wdV)

;;calculated weighted abundances (&errors) for each cluster within 2500 & 500

files ='a'
readcol, 'clu_names', files, format='A'
n_clu = n_elements(files)


for i=0,n_clu-1 do begin
    a_mw_2500 = -999 & a_em_2500 = -999 
    a_mw_500 = -999 & a_em_500 = -999
    a_mw_lim_2500 = -999 & a_em_lim_2500 = -999
    a_mw_lim_500 = -999 & a_em_lim_500 = -999

    restore, files[i]+'_basic.sav'
    x = where(abunderr eq 0.0,c)
    if (x[0] ne -1) then for j=0,c-1 do abunderr(x[j]) = abunderr(x[j]-1)
    restore, files[i]+'_derive.sav'

    if (r2500 gt 0) then begin
        x = where(rin_mpc lt r2500,c)
        rin_int = rin_mpc[x]
        rout_int = [rout_mpc[x(0:c-2)],r2500]
;mass-weighted

        int_array1 = gas_m_int(rin_int, rout_int, ne_best_fit)
        a_mw_2500 = total(abund[x]*int_array1)/total(int_array1)

        int_array2 = em_int(rin_int, rout_int, ne_best_fit)
        a_em_2500 = total(abund[x]*int_array2)/total(int_array2)

    ;;Now do error in both.  Assume error dominated by error in
    ;;abundance measurement
    
        rand = randomu(1000L,500,n_m)
        a_mw_dist_2500 = fltarr(500)
        a_em_dist_2500 = fltarr(500)
        for j=0,499 do begin
            abundnew = abund + (rand[j,*]*2.0 - 1.0)*abunderr
            a_mw_dist_2500[j] = total(abundnew[x]*int_array1)/$
              total(int_array1)
            a_em_dist_2500[j] = total(abundnew[x]*int_array2)/$
              total(int_array2)
        endfor
        a_mw_lim_2500 = get_limit(a_mw_2500, a_mw_dist_2500, 0.8413)
        a_em_lim_2500 = get_limit(a_em_2500, a_em_dist_2500, 0.8413)

    endif

    if (r500 gt 0) then begin
        x = where(rin_mpc lt r500,c)
        rin_int = rin_mpc[x]
        rout_int = [rout_mpc[x(0:c-2)],r500]
;mass-weighted

        int_array1 = gas_m_int(rin_int, rout_int, ne_best_fit)
        a_mw_500 = total(abund[x]*int_array1)/total(int_array1)

        int_array2 = em_int(rin_int, rout_int, ne_best_fit)
        a_em_500 = total(abund[x]*int_array2)/total(int_array2)

    ;;Now do error in both.  Assume error dominated by error in
    ;;abundance measurement
        
        rand = randomu(1000L,500,n_m)
        a_mw_dist_500 = fltarr(500)
        a_em_dist_500 = fltarr(500)
        for j=0,499 do begin
            abundnew = abund + (rand[j,*]*2.0 - 1.0)*abunderr
            a_mw_dist_500[j] = total(abundnew[x]*int_array1)/total(int_array1)
            a_em_dist_500[j] = total(abundnew[x]*int_array2)/total(int_array2)
        endfor
        a_mw_lim_500 = get_limit(a_mw_500, a_mw_dist_500, 0.8413)
        a_em_lim_500 = get_limit(a_em_500, a_em_dist_500, 0.8413)
        
    endif

    save, a_mw_2500, a_em_2500, a_mw_500, a_em_500, a_mw_lim_2500, $
      a_em_lim_2500, a_mw_lim_500, a_em_lim_500, $
      filename=strcompress(name,/remove_all)+'_abund.sav'

endfor

end


pro make_wt_figs
set_plot, 'ps'
device, filename='r500_2500_tabund.ps', /inches, ysize=9.0, yoffset=1.0


files ='a'
readcol, 'clu_500', files, format='A'
n_clu = n_elements(files)

amw_all = fltarr(n_clu,3)
aem_all = fltarr(n_clu,3)
tx_all = fltarr(n_clu,3)

for i=0,n_clu-1 do begin
    restore, 'part/'+files[i]+'_500derive.sav_part'
    tx_all[i,*] = [tx_500,tx_500_lim[0],tx_500_lim[1]]

    restore, files[i]+'_abund.sav'
    amw_all[i,*] = [a_mw_500, a_mw_lim_500[0], a_mw_lim_500[1]]
    aem_all[i,*] = [a_em_500, a_em_lim_500[0], a_em_lim_500[1]]
endfor

!p.multi=[0,1,2]
plot, tx_all[*,0], amw_all[*,0], psym=4, yrange=[0,1.0], xrange=[0,15], $
  ytitle='Abundance (mass-weighted)', title='R500 data'
for i=0,n_clu-1 do oplot, [tx_all[i,1],tx_all[i,2]], $
  [amw_all[i,0],amw_all[i,0]]
for i=0,n_clu-1 do oplot, [tx_all[i,0],tx_all[i,0]], $
  [amw_all[i,1],amw_all[i,2]]

plot, tx_all[*,0], aem_all[*,0], psym=4, yrange=[0,1.0], xrange=[0,15], $
  ytitle='Abundance (emission measure weighted)', xtitle='Temp (keV)'
for i=0,n_clu-1 do oplot, [tx_all[i,0],tx_all[i,0]], $
  [aem_all[i,1],aem_all[i,2]]
for i=0,n_clu-1 do oplot, [tx_all[i,1],tx_all[i,2]], $
  [aem_all[i,0],aem_all[i,0]]



readcol, 'clu_2500', files, format='A'
n_clu = n_elements(files)

amw_all = fltarr(n_clu,3)
aem_all = fltarr(n_clu,3)
tx_all = fltarr(n_clu,3)

for i=0,n_clu-1 do begin
    restore, 'part/'+files[i]+'_2500derive.sav_part'
    tx_all[i,*] = [tx_2500,tx_2500_lim[0],tx_2500_lim[1]]

    restore, files[i]+'_abund.sav'
    amw_all[i,*] = [a_mw_2500, a_mw_lim_2500[0], a_mw_lim_2500[1]]
    aem_all[i,*] = [a_em_2500, a_em_lim_2500[0], a_em_lim_2500[1]]
endfor

!p.multi=[0,1,2]
plot, tx_all[*,0], amw_all[*,0], psym=4, yrange=[0,1.0], xrange=[0,15], $
  ytitle='Abundance (mass-weighted)', title='R2500 data'
for i=0,n_clu-1 do oplot, [tx_all[i,1],tx_all[i,2]], $
  [amw_all[i,0],amw_all[i,0]]
for i=0,n_clu-1 do oplot, [tx_all[i,0],tx_all[i,0]], $
  [amw_all[i,1],amw_all[i,2]]

plot, tx_all[*,0], aem_all[*,0], psym=4, yrange=[0,1.0], xrange=[0,15], $
  ytitle='Abundance (emission measure weighted)', xtitle='Temp (keV)'
for i=0,n_clu-1 do oplot, [tx_all[i,0],tx_all[i,0]], $
  [aem_all[i,1],aem_all[i,2]]
for i=0,n_clu-1 do oplot, [tx_all[i,1],tx_all[i,2]], $
  [aem_all[i,0],aem_all[i,0]]


device, /close
set_plot, 'x'
end

    

pro make_avgzdist_fig
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
    restore, file_name[k]+'_basic.sav'
    x = where(abunderr eq 0.0,c)
    if (x[0] ne -1) then for j=0,c-1 do abunderr(x[j]) = abunderr(x[j]-1)
    restore, file_name[k]+'_derive.sav'

    for i=0,nbins-1 do begin
    ;type 1 -- fully contained
        x = where((rin_mpc/r500 ge rbin_low[i]) and $
                  (rout_mpc/r500 le rbin_hi[i]),c)
    
        if (c gt 0) then begin
            t_sum[i] = t_sum[i] + $
              total(abund[x]/abunderr[x]/abunderr[x])
            t_wt[i] = t_wt[i] + total(1.0/abunderr[x]/abunderr[x])
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
                  abund[x[j]]*wt/abunderr[x[j]]/abunderr[x[j]]
                t_wt[i] = t_wt[i] + wt/abunderr[x[j]]/abunderr[x[j]]
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
                  abund[x[j]]*wt/abunderr[x[j]]/abunderr[x[j]]
                t_wt[i] = t_wt[i] + wt/abunderr[x[j]]/abunderr[x[j]]
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
                  abund[x[j]]*wt/abunderr[x[j]]/abunderr[x[j]]
                t_wt[i] = t_wt[i] + wt/abunderr[x[j]]/abunderr[x[j]]
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
        if (c gt 0) then vals = [vals, abund[x]]

    ;type 2 -- in<binlow, out>binlow<binhi
        x = where((rin_mpc/r500 lt rbin_low[i]) and $
                  ((rout_mpc/r500 le rbin_hi[i]) and $
                   (rout_mpc/r500 gt rbin_low[i])),c)
        if (c gt 0) then vals = [vals, abund[x]]
       
    ;type 3 -- out>binout, in>binlow<binhi
        x = where((rout_mpc/r500 gt rbin_hi[i]) and $
                  ((rin_mpc/r500 lt rbin_hi[i]) and $
                   (rin_mpc/r500 ge rbin_low[i])),c)
        if (c gt 0) then vals = [vals, abund[x]]

    ;type 4 -- in<binlow, out>binhi
        x = where((rout_mpc/r500 gt rbin_hi[i]) and $
                  (rin_mpc/r500 lt rbin_low[i]),c)
        if (c gt 0) then vals = [vals, abund[x]]

    ;type 5 -- in&out<binlow - do nothing
    ;type 6 -- in&out>binhi - do nothing

    endfor
    
    vari[i] = sqrt(total((vals[1:*] - t_avg[i])^2)/(n_elements(vals)-2))
endfor


set_plot, 'ps'
device, filename='avgzdist_fig.ps'
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0 
!p.charthick=2.0

plot, rbin_mid, t_avg, psym=5, yrange=[0,0.8], xtitle='r/r!B500!N', $
  ytitle='Z (Z!B!9n!X!N)'
for j=0,nbins-1 do oplot, [rbin_mid[j],rbin_mid[j]], $
  t_avg[j]+[-t_avg_err[j],t_avg_err[j]]
for j=0,nbins-1 do oplot, [rbin_low[j],rbin_hi[j]], [t_avg[j],t_avg[j]]
oplot, rbin_mid, t_avg+vari, linestyle=1.0
oplot, rbin_mid, t_avg-vari, linestyle=1.0

device, /close
set_plot, 'x'

end
