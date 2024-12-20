    real ::        stepx,     stepy
    parameter ( stepx=0.1, stepy=0.1 ) !model x,y resolution
!      real        rhom
!      parameter ( rhom=1024. )  !mean density (critical density)
!      real        rtime
!      parameter ( rtime=0.000 ) !retension time (beaching hours)
    real ::        rin
    parameter ( rin=15000. )    !initial spreading radius (meters)
!      real        pi,                g
!      parameter ( pi = 3.1415926536, g = 9.81 )
    real ::        dg            !Degree ~ 1852m*60=111120m, 1853m*60=111180m
    parameter ( dg = 111317.0997 )

    real ::        kwind         !empirical factor for wind drift
    parameter ( kwind = 0.015) !from Yoon 2010

!     Kinematic viscosity
!     -------------------
    data vk / 1.e-6 /


!     particle position
    real :: xo,yo,zo
    real :: xn(nmax),yn(nmax),zn(nmax)

    real :: vk

    real :: xm,ym
                                             
    common/ibmxy/ xo(nmax),yo(nmax),zo(nmax)
    common/ibmxymax/xm,ym


