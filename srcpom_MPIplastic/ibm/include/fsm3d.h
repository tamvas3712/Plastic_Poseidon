!**********************************************************************
! fsm3d.h - last change: 10/04/97 (ank)
!        INCLUDED IN: fsm3d.f, output.f, sig2car.f, pom2grads.f
!**********************************************************************

    integer :: klev
    parameter (klev = 25)
    real :: dpt(klev)

    real :: fsm3d
    common / fsm3d_out / fsm3d(im,jm,klev)

!     Standard depths for output
!     --------------------------
    data dpt / 0., 5., 10., 15., 20., 25., 30., 35., 40., &
    50., 80., 100., 120., 160., 200., 240., 280., &
    340., 420., 500., 620., 850., 1250., 1750., 2000. /
