    PARAMETER (IM1=431,JM1=156) !original grid
    PARAMETER (IM=431,JM=156,KB=25) !modified grid to fit decomposition

    parameter(XLEFT=-7,YBOT=30.25,XRES=0.1,YRES=0.1)

    PARAMETER (KBM1=KB-1,IMM1=IM-1,JMM1=JM-1)

    real :: rbc
    integer :: lbc
    common / modelgrid / h(im,jm),alon(im,jm),alat(im,jm), &
    dx(im,jm),dy(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm), &
    z(kb),zz(kb),dz(kb),dzz(kb), &
    rbc(im,jm),lbc(im,jm),coast(im,jm)


