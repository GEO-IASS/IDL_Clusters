pro combine_r500
restore, 'mass_temp_res_corr.sav'
omega_m = 0.3 & omega_lam = 0.7
e_z =  sqrt(omega_m*(1.0 + zed)^3 + omega_lam)
masses = mass_new & j=2 

x = where((masses[*,1,0,0] gt 0.0) and $
          (strcompress(clu_name_mt, /remove_all) ne 'A2163'),nx)
y = where(masses[*,1,1,0] gt 0.0,ny)

;Here collect info for fit
fit_files = ['fit_res_orth','fit_res_yx','fit_res_fix']
fit_res = fltarr(4,2,3,2)

for i=0,2 do begin
    readcol, fit_files[i], b, m
    fit_res[*,0,i,0] = b[0:3]
    fit_res[*,0,i,1] = m[0:3]
    fit_res[*,1,i,0] = b[4:7]
    fit_res[*,1,i,1] = m[4:7]
endfor
color = [50,200,130]

set_plot, 'ps'
device, filename='combine_r500.ps', /inches, xsize=7.5, /color, ysize=6.8
!p.multi=[0,2,2]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

;Mass vs Temp r500
yaxis = e_z(y)*10.^(masses[y,1,1,0])
xaxis = temps[y,2,1,0]
xup = temps[y,2,1,2] & xdown = temps[y,2,1,1]
yup = e_z(y)*10.^(masses[y,1,1,2]) & ydown =e_z(y)*10.^(masses[y,1,1,1])

     plot, xaxis, yaxis, psym=4, xtitle='T!DX!N (keV)', /xlog, /ylog, $
       ytitle='E(z)!9X!XM!D500!N (M!D!9n!X!N)', xrange=[1.0,20.0], xstyle=1, $
       yrange=[5d13,5d15], ystyle=1, xmargin=[8,2], ymargin=[3,2]
     for i=0,ny-1 do oplot, [xaxis[i],xaxis[i]], [ydown[i],yup[i]]
     for i=0,ny-1 do oplot, [xdown[i],xup[i]], [yaxis[i],yaxis[i]]

     xfit = findgen(21)/20. * (10^!x.crange[1] - 10^!x.crange[0]) + $
       10^!x.crange[0]
     for i=0,2 do begin
         yfit = 10.0^fit_res[0,1,i,0] * (xfit/5.0)^fit_res[0,1,i,1]
         oplot, xfit, yfit, thick=5.0, color=color[i]
         ;legend, ['Orth','Y|X','Fix'], linestyle=[0,0,0],colors=color
     endfor

;Mass vs Y_X r500
yaxis = e_z(y)^(2./5.)*10.^(masses[y,1,1,0])
xaxis = temps[y,2,1,0]*10.^(masses[y,0,1,0])
xup = temps[y,2,1,2]*10.^(masses[y,0,1,2])
xdown = temps[y,2,1,1]*10.^(masses[y,0,1,1])
yup = e_z(y)^(2./5.)*10.^(masses[y,1,1,2]) 
ydown = e_z(y)^(2./5.)*10.^(masses[y,1,1,1])

     plot, xaxis, yaxis, psym=4, /xlog, /ylog, ymargin=[3,2], $
       xtitle='Y!DX!N !9:!X T!DX!NM!Dg,500!N (keV M!D!9n!X!N)', $
       ytitle='E(z)!U2/5!N!9X!XM!D500!N (M!D!9n!X!N)', xmargin=[7,3], $
       xrange=[2d13,1d16], yrange=[5d13,5d15], xstyle=1, ystyle=1
     for i=0,ny-1 do oplot, [xaxis[i],xaxis[i]], [ydown[i],yup[i]]
     for i=0,ny-1 do oplot, [xdown[i],xup[i]], [yaxis[i],yaxis[i]]

     xfit = findgen(21)/20. * (10^!x.crange[1] - 10^!x.crange[0]) + $
       10^!x.crange[0]
     for i=0,2 do begin
         yfit = 10.0^fit_res[1,1,i,0] * (xfit/3d14)^fit_res[1,1,i,1]
         oplot, xfit, yfit, thick=5.0, color=color[i]
         ;legend, ['Orth','Y|X','Fix'], linestyle=[0,0,0],colors=color
     endfor

;Gas vs Total Mass r500
yaxis = 10.^masses[y,1,1,0] & xaxis = 10.^masses[y,0,1,0]
xup = 10.^masses[y,0,1,2] & xdown = 10.^masses[y,0,1,1]
yup = 10.^masses[y,1,1,2] & ydown = 10.^masses[y,1,1,1]

     plot, xaxis, yaxis, psym=4, xtitle='M!Dg,500!N (M!D!9n!X!N)', $
       ytitle = 'M!D500!N (M!D!9n!X!N)', /xlog, /ylog, xmargin=[8,2], $
       xrange=[1d13,7d14], xstyle=1, yrange=[7d13,5d15], ystyle=1, $
       ymargin=[4,1]
     for i=0,ny-1 do oplot, [xaxis[i],xaxis[i]], [ydown[i],yup[i]]
     for i=0,ny-1 do oplot, [xdown[i],xup[i]], [yaxis[i],yaxis[i]]

     xfit = findgen(21)/20. * (10^!x.crange[1] - 10^!x.crange[0]) + $
       10^!x.crange[0]
     for i=0,2 do begin
         yfit = 10.0^fit_res[2,1,i,0] * (xfit/1d13)^fit_res[2,1,i,1]
         oplot, xfit, yfit, thick=5.0, color=color[i]
         ;legend, ['Orth','Y|X','Fix'], linestyle=[0,0,0],colors=color
     endfor

;Gas Fraction vs Total Mass r500
yaxis = 10.^(masses[y,0,1,0] - masses[y,1,1,0]) & xaxis = 10.^(masses[y,1,1,0])
xup = 10.^(masses[y,1,1,2]) & xdown = 10.^(masses[y,1,1,1])
yup = 10.^(masses[y,0,1,2] - masses[y,1,1,1])
ydown = 10.^(masses[y,0,1,1] - masses[y,1,1,2])

     plot, xaxis, yaxis, psym=4, xtitle='M!D500!N (M!D!9n!X!N)', $
       ytitle = 'f!Dg,500!N', /xlog, xmargin=[7,3], ymargin=[4,1], $
       xrange=[7d13,5d15], xstyle=1, yrange=[0,0.25], ystyle=1
     for i=0,ny-1 do oplot, [xaxis[i],xaxis[i]], [ydown[i],yup[i]]
     for i=0,ny-1 do oplot, [xdown[i],xup[i]], [yaxis[i],yaxis[i]]

     xfit = findgen(21)/20. * (!x.crange[1] - !x.crange[0]) + !x.crange[0]
     i=0
     yfit = fit_res[3,1,i,0] + fit_res[3,1,i,1]*(xfit - 15.0)
     oplot, 10.^xfit, yfit, thick=5.0, color=color[i]
     ;legend, ['Orth','Y|X'], linestyle=[0,0],colors=color[0:1]

device, /close
set_plot, 'x'

end



pro plot_scat_pap
restore, 'scat_res_2.sav'

x = where(masses[*,1,1,0] gt 0.0,nx)
r = 1
label = 'r500'


set_plot,'ps'


;First do for M-T
device, filename='plot_scat_mt_v2.ps', /inches, xsize=7.5, /color, ysize=6.8
!p.multi=[0,2,2]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

yaxis = scat_m_t[x,1,r]
yerr_d = 10.^(masses[x,1,r,1] - masses[x,1,r,0]) - 1.0
yerr_u = 10.^(masses[x,1,r,2] - masses[x,1,r,0]) - 1.0
yrange = [min(yaxis-yerr_d),max(yaxis+yerr_u)]

;vs ellip
xaxis = ellip[x,0] & xerr = ellip[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Ellipticity', $
       ytitle='!4d!XM(T!DX!N)/M', yrange=yrange, xrange=[0,0.15], $
       xmargin=[8,2], ymargin=[3,2]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

;vs Asym 
xaxis = asym[x,0] & xerr = asym[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Asymmetry', $
       ytitle='!4d!XM(T!DX!N)/M', yrange=yrange, xmargin=[7,3], ymargin=[3,2]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

;vs beta
xaxis = beta_t_res[x,r,0]+beta_ne_res[x,r,0]
xup = beta_t_res[x,r,2]+beta_ne_res[x,r,2]
xdown = beta_t_res[x,r,1]+beta_ne_res[x,r,1]

     plot, xaxis, yaxis, psym=4, xtitle='!4b!X!Deff!N+!4b!X!Dt!N', $
       ytitle='!4d!XM(T!DX!N)/M', yrange=yrange, xrange=[0,3.0], $
       xmargin=[8,2], ymargin=[4,1]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

;vs conc
xaxis = conc[x,r,0] & xdown = conc[x,r,1] & xup = conc[x,r,2]

     plot, xaxis, yaxis, psym=4, xtitle='Concentration', $
       ytitle='!4d!XM(T!DX!N)/M', yrange=yrange, xrange=[0,25], $
       xmargin=[7,3], ymargin=[4,1]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

device, /close

device, filename='plot_scat_my.ps', /inches, xsize=7.5, /color, ysize=6.8
!p.multi=[0,2,2]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

;M-Y
yaxis = scat_m_y[x,1,r]
yerr_d = 10.^(masses[x,1,r,1] - masses[x,1,r,0]) - 1.0
yerr_u = 10.^(masses[x,1,r,2] - masses[x,1,r,0]) - 1.0
yrange = [min(yaxis-yerr_d),max(yaxis+yerr_u)]

;vs beta
xaxis = beta_t_res[x,r,0]+beta_ne_res[x,r,0]
xup = beta_t_res[x,r,2]+beta_ne_res[x,r,2]
xdown = beta_t_res[x,r,1]+beta_ne_res[x,r,1]

     plot, xaxis, yaxis, psym=4, xtitle='!4b!X!Deff!N+!4b!X!Dt!N', $
       ytitle='!4d!XM(Y!DX!N)/M', yrange=yrange, xrange=[0,3.0], $
       xmargin=[8,2], ymargin=[3,2]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

;vs conc
xaxis = conc[x,r,0] & xdown = conc[x,r,1] & xup = conc[x,r,2]

     plot, xaxis, yaxis, psym=4, xtitle='Concentration', $
       ytitle='!4d!XM(Y!DX!N)/M', yrange=yrange, xrange=[0,25], $
       xmargin=[7,3], ymargin=[3,2]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

;vs ellip
xaxis = ellip[x,0] & xerr = ellip[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Ellipticity', $
       ytitle='!4d!XM(Y!DX!N)/M', yrange=yrange, xrange=[0,0.15], $
       xmargin=[8,2], ymargin=[4,1]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

;vs Asym 
xaxis = asym[x,0] & xerr = asym[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Asymmetry', $
       ytitle='!4d!XM(Y!DX!N)/M', yrange=yrange, xmargin=[7,3], ymargin=[4,1]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

device, /close

device, filename='plot_scat_mmg.ps', /inches, xsize=7.5, /color, ysize=6.8
!p.multi=[0,2,2]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

;M-Mg
yaxis = scat_m_mg[x,1,r]
yerr_d = 10.^(masses[x,1,r,1] - masses[x,1,r,0]) - 1.0
yerr_u = 10.^(masses[x,1,r,2] - masses[x,1,r,0]) - 1.0
yrange = [min(yaxis-yerr_d),max(yaxis+yerr_u)]

;vs beta
xaxis = beta_t_res[x,r,0]+beta_ne_res[x,r,0]
xup = beta_t_res[x,r,2]+beta_ne_res[x,r,2]
xdown = beta_t_res[x,r,1]+beta_ne_res[x,r,1]

     plot, xaxis, yaxis, psym=4, xtitle='!4b!X!Deff!N+!4b!X!Dt!N', $
       ytitle='!4d!XM(M!Dg!N)/M', yrange=yrange, xrange=[0,3.0], $
       xmargin=[8,2], ymargin=[3,2]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

;vs conc
xaxis = conc[x,r,0] & xdown = conc[x,r,1] & xup = conc[x,r,2]

     plot, xaxis, yaxis, psym=4, xtitle='Concentration', $
       ytitle='!4d!XM(M!Dg!N)/M', yrange=yrange, xrange=[0,25], $
       xmargin=[7,3], ymargin=[3,2]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

;vs ellip
xaxis = ellip[x,0] & xerr = ellip[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Ellipticity', $
       ytitle='!4d!XM(M!Dg!N)/M', yrange=yrange, xrange=[0,0.15], $
       xmargin=[8,2], ymargin=[4,1]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

;vs Asym 
xaxis = asym[x,0] & xerr = asym[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Asymmetry', $
       ytitle='!4d!XM(M!Dg!N)/M', yrange=yrange, xmargin=[7,3], ymargin=[4,1]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

device, /close
set_plot, 'x'
end


pro beta_figure
restore, 'scat_res_2.sav'

x = where(masses[*,1,1,0] gt 0.0,nx)
r = 1
label = 'r500'


set_plot,'ps'
device, filename='beta_figure.ps', /inches, xsize=3.7, ysize=4.5
!p.multi=[0,1,2]
!p.charsize=1.4 & !p.thick=2.1 & !x.thick=2.1 & !y.thick=2.1
!p.charthick=2.5

plot, temps[x,2,r,0], beta_t_res[x,r,0], psym=5, ytitle='!4b!X!Dt!N', $
  xtickname=[' ',' ',' ',' ',' ',' ',' '], ymargin=[-1,2], yrange=[-0.4,2.0], $
  ystyle=1, xrange=[1.0,15], xstyle=1, xmargin=[8,3]
for i=0,nx-1 do oplot, temps[x(i),2,r,0]+[0,0], $
  [beta_t_res[x(i),r,1],beta_t_res[x(i),r,2]]
oplot, [0,20], [0.17,0.17], linestyle=1, thick=5, color=130
oplot, [0,20], [mean(beta_t_res[x,r,0]),mean(beta_t_res[x,r,0])], thick=5.0, $
  color=130

plot, temps[x,2,r,0], beta_ne_res[x,r,0], psym=6, ytitle='!4b!X!Deff!N', $
  yrange=[0.2,0.99], ystyle=1, xrange=[1.0,15], xstyle=1, ymargin=[4,1], $
  xmargin=[8,3], xtitle='T!DX!N (keV)'
for i=0,nx-1 do oplot, temps[x(i),2,r,0]+[0,0], $
  [beta_ne_res[x(i),r,1],beta_ne_res[x(i),r,2]]
oplot, [0,20], [0.78,0.78], linestyle=1, thick=5, color=130
oplot, [0,20], [mean(beta_ne_res[x,r,0]),mean(beta_ne_res[x,r,0])], $
  thick=5.0, color=130


device, /close
set_plot,'x'
end


pro temp_dist_fig
restore, 'scat_res_2.sav'

x = where(masses[*,1,1,0] gt 0.0,nx)
beta_facs = moment(beta_t_res[x,1,0])
sig = 1.0

color = [50,130,200]

set_plot,'ps'
device, filename='temp_dists.ps', /color
!p.multi=[0,1,1]
loadct,6 

;Do first for i=0
restore, strcompress(clu_name_mt[x(0)],/remove_all)+'_vars.sav'
r_all = [r_mpc,max(rout_mpc)]/radii_res[1]
t = temp_dist(r_all*radii_res[1], temp_best_fit)

plot, r_all, t/temp_res[2,1], /xlog, yrange=[0,2.0]

if (beta_t_res[x(0),1,0] ge beta_facs[0]-sig*sqrt(beta_facs[1])) and $
  (beta_t_res[x(0),1,0] le beta_facs[0]+sig*sqrt(beta_facs[1])) then $
  col = color[0] $
else if (beta_t_res[x(0),1,0] lt beta_facs[0]-sig*sqrt(beta_facs[1])) $
  then col = color[1] $
else if (beta_t_res[x(0),1,0] gt beta_facs[0]+sig*sqrt(beta_facs[1])) $
  then col = color[2]

oplot, r_all, t/temp_res[2,1], color=col

for i=1,nx-1 do begin    
    restore, strcompress(clu_name_mt[x(i)],/remove_all)+'_vars.sav'
    r_all = [r_mpc,max(rout_mpc)]/radii_res[1]
    t = temp_dist(r_all*radii_res[1], temp_best_fit)

if (beta_t_res[x(i),1,0] ge beta_facs[0]-sig*sqrt(beta_facs[1])) and $
  (beta_t_res[x(i),1,0] le beta_facs[0]+sig*sqrt(beta_facs[1])) then $
  col = color[0] $
else if (beta_t_res[x(i),1,0] lt beta_facs[0]-sig*sqrt(beta_facs[1])) $
  then col = color[1] $
else if (beta_t_res[x(i),1,0] gt beta_facs[0]+sig*sqrt(beta_facs[1])) $
  then col = color[2]

    oplot, r_all, t/temp_res[2,1], color=col
    
endfor

device, /close
set_plot,'x'

end


;New scatter plots
pro plot_scats
restore, 'scat_res_2.sav'

x = where(masses[*,1,1,0] gt 0.0,nx)
r = 1
label = 'r500'


set_plot,'ps'

;Just do scatter for Yx ellip/asym as example.
device, filename='plot_yx_e+a.ps', /inches, xsize=4.0, /color, ysize=6.5
!p.multi=[0,1,2]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

yaxis = scat_m_y[x,1,r]
yerr_d = 10.^(masses[x,1,r,1] - masses[x,1,r,0]) - 1.0
yerr_u = 10.^(masses[x,1,r,2] - masses[x,1,r,0]) - 1.0
yrange = [min(yaxis-yerr_d),max(yaxis+yerr_u)]

;vs ellip
xaxis = ellip[x,0] & xerr = ellip[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Ellipticity', $
       ytitle='!4d!XM(Y!DX!N)/M', yrange=yrange, xrange=[0,0.15], $
       xmargin=[9,3], ymargin=[3,2]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

;vs Asym 
xaxis = asym[x,0] & xerr = asym[x,1]

     plot, xaxis, yaxis, psym=4, xtitle='Asymmetry', $
       ytitle='!4d!XM(Y!DX!N)/M', yrange=yrange, xmargin=[9,3], ymargin=[4,1]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, xaxis[j]+[-xerr[j],xerr[j]], [yaxis[j],yaxis[j]]

device, /close


;Also do scatter for T, Mg, Yx vs beta_t + beta_eff
device, filename='plot_beta_y+mg+t.ps', /inches, xsize=4.0, /color, ysize=7.0
!p.multi=[0,1,3]
!p.charsize=2.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0
loadct, 6

;vs beta
xaxis = beta_t_res[x,r,0]+beta_ne_res[x,r,0]
xup = beta_t_res[x,r,2]+beta_ne_res[x,r,2]
xdown = beta_t_res[x,r,1]+beta_ne_res[x,r,1]

;M-T
yaxis = scat_m_t[x,1,r]
yerr_d = 10.^(masses[x,1,r,1] - masses[x,1,r,0]) - 1.0
yerr_u = 10.^(masses[x,1,r,2] - masses[x,1,r,0]) - 1.0

     plot, xaxis, yaxis, psym=4, ytitle='!4d!XM(T!DX!N)/M', $
       yrange=[-1.6,0.8], ystyle=1, $
       xrange=[0,3.0], ymargin=[0,2], xtickname=[' ',' ',' ',' ',' ']
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

;M-Mg
yaxis = scat_m_mg[x,1,r]
yerr_d = 10.^(masses[x,1,r,1] - masses[x,1,r,0]) - 1.0
yerr_u = 10.^(masses[x,1,r,2] - masses[x,1,r,0]) - 1.0
yrange = [min(yaxis-yerr_d),max(yaxis+yerr_u)]

     plot, xaxis, yaxis, psym=4, ytitle='!4d!XM(M!Dg!N)/M', $
       yrange=[-0.39,0.58], ystyle=1, $
       xrange=[0,3.0], xtickname=[' ',' ',' ',' ',' '], ymargin=[2,0]
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

;M-Yx
yaxis = scat_m_y[x,1,r]
yerr_d = 10.^(masses[x,1,r,1] - masses[x,1,r,0]) - 1.0
yerr_u = 10.^(masses[x,1,r,2] - masses[x,1,r,0]) - 1.0
yrange = [min(yaxis-yerr_d),max(yaxis+yerr_u)]

     plot, xaxis, yaxis, psym=4, xtitle='!4b!X!Deff!N+!4b!X!Dt!N', $
       ytitle='!4d!XM(Y!DX!N)/M', yrange=[-0.7,0.7], xrange=[0,3.0], $
       ymargin=[4,-2], ystyle=1
     for j=0,nx-1 do oplot, [xaxis[j],xaxis[j]], yaxis[j]+[yerr_d[j],yerr_u[j]]
     for j=0,nx-1 do oplot, [xup[j],xdown[j]], [yaxis[j],yaxis[j]]

device, /close
set_plot, 'x'
end


pro fg_dist
restore, 'mass_temp_res_corr.sav'
masses = mass_new 

y = where(masses[*,1,1,0] gt 0.0,ny)

fg = 10.^(masses[y,0,1,0] - masses[y,1,1,0]) 

nbins=25

hist = histogram(fg, min=0.0, max=0.25, nbins=nbins)
fg_bin = (findgen(nbins)+0.5)/(nbins-1) * 0.25


set_plot,'ps'
device, filename='plot_fg_dist.ps', /inches, xsize=4.0, /color, ysize=3.5
!p.multi=[0,1,1]
!p.charsize=1.3 & !p.thick=2.0 & !x.thick=2.0 & !y.thick=2.0
!p.charthick=3.0 & !x.ticklen=0.04 & !y.ticklen=0.03
loadct, 6

plot, fg_bin, hist, psym=10, xrange=[0.05,0.20], xstyle=1, yrange=[0,10], $
  xtitle = 'f!Dg,500!N', ytitle='Number', xmargin=[7,3], xticks=3, xminor=5

device, /close
set_plot, 'x'

end




