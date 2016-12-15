pro plot_parms_a,nm, parm

;The parameters are (as given by Vikhlinin) .
;    0 = T_0
;    1 = r_cool
;    2 = a_cool
;    3 = T_min
;    4 = r_t
;    5 = a
;    6 = b
;    7 = c
!x.range=[1,15]
!y.range=[-0.2,0.2]
!x.style=1
!y.style=1
plot,parm[*,0],parm[*,5],psym=1, xtitle="T0(KeV)",ytitle='fitted a'
;oplot,parm[*,0],parm[*,6],psym=2
;oplot,parm[*,0],parm[*,7],psym=4



return
end
