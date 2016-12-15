pro plot_parms,nm, parm

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
!x.style=1
!y.style=1
!p.multi=[0,2,3,0,0]
!y.range=[-1.,6.]
plot,parm[*,0],parm[*,7]/parm[*,6],psym=1, xtitle="T0(KeV)",ytitle='fitted c/b'
!y.range=[-1,3.]
plot,parm[*,0],-1*parm[*,5]/parm[*,7],psym=1, xtitle="T0(KeV)",ytitle='asymptote (-a/c)'
!y.range=[-0.2,0.2]
plot,parm[*,0],parm[*,5],psym=1, xtitle="T0(KeV)",ytitle='fitted a'
!y.range=[-2,6]
plot,parm[*,0],parm[*,6],psym=1, xtitle="T0(KeV)",ytitle='fitted b'
!y.range=[-2,6]
plot,parm[*,0],parm[*,7],psym=1, xtitle="T0(KeV)",ytitle='fitted c'
;oplot,parm[*,0],parm[*,6],psym=2
;oplot,parm[*,0],parm[*,7],psym=4



return
end
