function get_limit, best, sample, per
;upper limit
x = where(sample ge best)
if (x[0] eq -1) then ulim = -999 else begin
    a = sort(sample[x])
    spot = per*n_elements(x)
    ulim = sample[x(a[spot])]
endelse

;lower limit
x = where(sample le best)
if (x[0] eq -1) then llim = -999 else begin
    a = reverse(sort(sample[x]))
    spot = per*n_elements(x)
    llim = sample[x(a[spot])]
endelse

lim = [llim,ulim]
return, lim
end


function temp_dist, X, P

z = (X/P[1])^(P[2])
tc = (z + P[3])/(1.d + z)

w = (X/P[4])^(-1.d*P[5])
v = (X/P[4])^(P[6])
to = w/(1.d + v)^(P[7]/P[6])

tmod = P[0]*tc*to

return, tmod
end


function ne_dist, r, a

cusp = (r/a[1])^(-a[3])
brk = 1.0/(1.0 + r^3/a[4]^3)^(a[5]/3.) 
king = 1.0/(1.0 + r^2/a[1]^2)^(3*a[2] - a[3]/2.)
k2 = 1.0/(1.0 + r^2/a[7]^2)^(3*a[8])

ne2_fxn = a[0]*a[0]*cusp*brk*king + a[6]*a[6]*k2

ne_fxn = sqrt(ne2_fxn)

return, ne_fxn
end


function rhog_dist, r, a
;Do in solar mass/mpc^3
mp = 1.67262158d-27; kg
msol = 1.98892d30 ;kg
mpc_2_cm = 3.086d24

ne_d = ne_dist(r,a)

rho_g = sqrt(1.624)*mp*ne_d/msol*(mpc_2_cm*mpc_2_cm*mpc_2_cm)

return, rho_g
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


function gas_mass, r, a
common int_fun, par

par = a
mp = 1.67262158d-27; kg
msol = 1.98892d30 ;kg
mpc_2_cm = 3.086d24

nr = n_elements(r)
int = fltarr(nr)
for i=0, nr-1 do int[i] = qromb('int_ner2', 1e-4, r[i], jmax=25)

log_m_gas = alog10(4*!pi*sqrt(1.624)*mp) + alog10(int) + 3*alog10(mpc_2_cm) $
  - alog10(msol)

return, log_m_gas
end


function total_mass, r, ne_par, temp_par

t = temp_dist(r, temp_par)

dnpne_dr = -6.0*ne_par[6]^2*ne_par[8]/ne_par[7]^2 * $
  r/(1.0 + r^2/ne_par[7]^2)^(3*ne_par[8] + 1.0) - ne_par[0]^2 * $
  (r/ne_par[1])^(-ne_par[3]) * 1.0/(1.0 + r^3/ne_par[4]^3)^(ne_par[5]/3.) * $
  1.0/(1.0 + r^2/ne_par[1]^2)^(3*ne_par[2] - ne_par[3]/2.) * $
  (ne_par[3]/r + 2.0*r*(3*ne_par[2] - ne_par[3]/2.)/(ne_par[1]^2 * $
  (1.0+ r^2/ne_par[1]^2)) + ne_par[5]*r^2/(ne_par[4]^3 * $
  (1.0 + r^3/ne_par[4]^3)))

sqrt_npne = ne_dist(r,ne_par)
dlogp_dlogr = r/(2.0*sqrt_npne^2) * dnpne_dr

dlogt_dlogr = temp_par[2]*(r/temp_par[1])^(temp_par[2])*(1.0 - temp_par[3])/$
  ((1.d + (r/temp_par[1])^(temp_par[2]))*((r/temp_par[1])^(temp_par[2]) + $
  temp_par[3])) - temp_par[5] - temp_par[7]*(r/temp_par[4])^(temp_par[6])/$
  (1.d +(r/temp_par[4])^(temp_par[6]))

tot_mass = alog10(-3.68d13*t*r*(dlogp_dlogr + dlogt_dlogr)) 

return, tot_mass
end


function entropy_dist, r, ne_par, temp_par
t = temp_dist(r, temp_par)

ne_d = ne_dist(r, ne_par)

ent = t/ne_d^(2./3.)

return, ent
end


function mass_comp, r, nep, tep, del, z
msol = 1.98892d30 ;kg
mpc_2_m = 3.086d22
G = 6.673d-11 ;m^3 /kg/s^2
omega_m = 0.3
omega_lam = 0.7
h0 = 70.0 ;km/s/Mpc
ez2 = omega_m*(1.0 + z)^3 + omega_lam
hz2 = h0*h0*ez2

mass_real = total_mass(r, nep, tep)

mass_cr = alog10(hz2*del*r^3/(2*G) *1000.0^2) + alog10(mpc_2_m) - alog10(msol)

return, mass_real - mass_cr
end


function calc_temp_spec_v2, mode, calfile, rin, rout
common tspec_pars, tcal, fcont_cal, fline_cal, emean_cal, x1, x2, ncal, acont

;Set parameters for Chandra by default
acont = 0.875 & beta = 1 & delta1 = 0.19 & delta2 = 0.25
;If instrument setup is different then change
if (mode eq 'xmm/pn') then begin
    acont = 0.79 & beta = 0.75 & delta1 = 0.270 & delta2 = 0.225
endif 
if  (mode eq 'xmm/mos') then begin
    acont = 0.90 & beta = 1 & delta1 = 0.19 & delta2 = 0.22
endif
if (mode eq 'xmm/mos+pn') then begin
    acont = 0.91 & beta = 0.90 & delta1 = 0.19 & delta2 = 0.21
endif
if (mode eq 'asca/sis') then begin
    acont = 0.875 & beta = 0.80 & delta1 = 0.20 & delta2 = 0.22
endif
if (mode eq 'asca/gis') then begin
    acont = 0.79 & beta = 0.75 & delta1 = 0.26 & delta2 = 0.30
endif

;read in calfile
readcol, calfile, tcal, fcont_cal, fline_cal, emean_cal, x1, x2, /silent
ncal = n_elements(tcal)

fluxline = int_2d('ltzp2dv',[rin, rout], 'z_lim', 96)
emean_n = int_2d('ftltzp2dv',[rin, rout], 'z_lim', 96)/fluxline
Tline_n = interpol(tcal, emean_cal, emean_n)

Tcont = int_2d('ctp2tatdv',[rin, rout], 'z_lim', 96)/$
  int_2d('ctp2tadv',[rin, rout], 'z_lim', 96)
fluxcont = int_2d('ctp2dv',[rin, rout], 'z_lim', 96)

f_cont = fluxcont/(fluxcont + fluxline)
f_line = 1 - f_cont
x = exp(-1*(f_line/delta1)^(2.0*beta)) * exp(-1*(f_line/delta2)^8.0)

temp_spec = x*Tcont + (1.0 - x)*Tline_n

return, temp_spec
end

function ltzp2dv, p, z
common int_issues, rmax, ne_pars, temp_pars, abund, rin_mpc, rout_mpc
common tspec_pars, tcal, fcont_cal, fline_cal, emean_cal, x1, x2, ncal, alpha
r = sqrt(p^2 + z^2)
t = temp_dist(r, temp_pars)
sqrt_nenp = ne_dist(r, ne_pars)

l_t = interpol(fline_cal, tcal, t)

n = n_elements(r) & ab = fltarr(n) & n_m = n_elements(abund)
for i=0,n-1 do begin
    x = where((r[i] ge rin_mpc) and (r[i] le rout_mpc))
    if (x[0] eq -1) then $
      if (r[i] gt rout_mpc[n_m-1]) then ab[i] = abund[n_m-1] else stop $
    else ab[i] = abund(x)
endfor

return, p * l_t*ab*sqrt_nenp*sqrt_nenp
end

function ftltzp2dv, p, z
common int_issues, rmax, ne_pars, temp_pars, abund, rin_mpc, rout_mpc
common tspec_pars, tcal, fcont_cal, fline_cal, emean_cal, x1, x2, ncal, alpha
r = sqrt(p^2 + z^2)
t = temp_dist(r, temp_pars)
sqrt_nenp = ne_dist(r, ne_pars)

l_t = interpol(fline_cal, tcal, t)
emean = interpol(emean_cal, tcal, t)

n = n_elements(r) & ab = fltarr(n) & n_m = n_elements(abund)
for i=0,n-1 do begin
    x = where((r[i] ge rin_mpc) and (r[i] le rout_mpc))
    if (x[0] eq -1) then $
      if (r[i] gt rout_mpc[n_m-1]) then ab[i] = abund[n_m-1] else stop $
    else ab[i] = abund(x)
endfor

return, p * emean*l_t*ab*sqrt_nenp*sqrt_nenp
end

function ctp2tatdv, p, z
common int_issues, rmax, ne_pars, temp_pars
common tspec_pars, tcal, fcont_cal, fline_cal, emean_cal, x1, x2, ncal, alpha

r = sqrt(p^2 + z^2)
t = temp_dist(r, temp_pars)
ct = interpol(fcont_cal, tcal, t)

sqrt_nenp = ne_dist(r, ne_pars)

return, p * ct*sqrt_nenp*sqrt_nenp*t^(1.0 - alpha)
end

function ctp2tadv, p, z
common int_issues, rmax, ne_pars, temp_pars
common tspec_pars, tcal, fcont_cal, fline_cal, emean_cal, x1, x2, ncal, alpha

r = sqrt(p^2 + z^2)
t = temp_dist(r, temp_pars)
ct = interpol(fcont_cal, tcal, t)

sqrt_nenp = ne_dist(r, ne_pars)

return, p * ct*sqrt_nenp*sqrt_nenp*t^(-1*alpha)
end

function ctp2dv, p, z
common int_issues, rmax, ne_pars, temp_pars
common tspec_pars, tcal, fcont_cal, fline_cal, emean_cal, x1, x2, ncal, alpha

r = sqrt(p^2 + z^2)
t = temp_dist(r, temp_pars)
ct = interpol(fcont_cal, tcal, t)

sqrt_nenp = ne_dist(r, ne_pars)

return, p * ct*sqrt_nenp*sqrt_nenp
end

function z_lim, p
common int_issues, rmax, ne_pars, temp_pars
return, [0.0, rmax]
end
