!    ECOLOGY....
!**********************************************************************
!    ecooutput.h
!               INCLUDED IN: output.f
!**********************************************************************

    real :: n1pmm(N_COMP), n3nmm(N_COMP), n4nmm(N_COMP), n5smm(N_COMP), &
    p1cmm(N_COMP), p2cmm(N_COMP), p3cmm(N_COMP), z5cmm(N_COMP), &
    z6cmm(N_COMP), b1cmm(N_COMP), r6cmm(N_COMP), r1cmm(N_COMP), &
    p4cmm(N_COMP), z4cmm(N_COMP), chlmm(N_COMP), &
    o2omm(N_COMP), netppmm(N_COMP), prodB1mm(N_COMP), &
    xepsmm(N_COMP), eirmm(N_COMP),o3cmm(N_COMP),o4nmm(N_COMP), &
    totamm(N_COMP) ,phxmm(N_COMP), carbamm(N_COMP),bicarbmm(N_COMP), &
    Carbmm(N_COMP),pco2wmm(N_COMP),cco2mm(N_COMP),fairmm(N_COMP)
    real :: q6pmm(I_COMPb:N_COMP),q1pmm(I_COMPb:N_COMP), &
    q6cmm(I_COMPb:N_COMP),q6nmm(I_COMPb:N_COMP), &
    q6smm(I_COMPb:N_COMP),g2omm(I_COMPb:N_COMP)

    real :: elemm(im,jm),tmm(im,jm,kb),smm(im,jm,kb),akhmm(im,jm,kb), &
    umm(im,jm,kb),vmm(im,jm,kb)

    common / eco_outarr / n1pmm, n3nmm, n4nmm, n5smm, p1cmm, p2cmm, &
    p3cmm, z5cmm, z6cmm, b1cmm, r6cmm, r1cmm, &
    p4cmm, z4cmm, chlmm, o2omm, netppmm,o3cmm,o4nmm, &
    prodB1mm, xepsmm, eirmm,q6pmm,q1pmm, &
    elemm, tmm, smm, akhmm, umm, vmm, &
    totamm, phxmm, carbamm, bicarbmm,Carbmm, pco2wmm, cco2mm,fairmm
                          
