    parameter(im2=872,jm2=317,kb2=25, &
    im=431,jm=156,kb=25, &
    kmat=5,kl=105)

    COMMON/FINE/h(im,jm),alon(im,jm),alat(im,jm), &
    dx(im,jm),dy(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm), &
    z(kb),zz(kb),dz(kb),dzz(kb)

    COMMON/COARSE/h2(im2,jm2),alon2(im2,jm2),alat2(im2,jm2), &
    dx2(im2,jm2),dy2(im2,jm2),fsm2(im2,jm2),dum2(im2,jm2), &
    dvm2(im2,jm2), &
    z2(kb2),zz2(kb2),dz2(kb2),dzz2(kb2)


