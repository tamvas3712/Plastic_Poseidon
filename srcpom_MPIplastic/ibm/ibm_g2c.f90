# 1 "ibm_g2c.F90"
subroutine Geo2Cart (xc, yc, slon, slat, xg, yg, pi, dg)
!! Converts Geographic coordinates to Cartesian (in meters)
    real :: xc, yc          !Output fields: Cartesian coordinates
    real :: slon, slat      !Starting Lon/Lat grid coordinates
    real :: xg, yg          !Geographic coordinates
    real :: aa, pi, dg      !work field aa, pi=3.14..., dg=1degree in meters

    yc = (yg - slat) * dg
    aa = pi * yg / 180.
    xc = (xg - slon) * (dg*cos(aa))

!      print *,'Cartesian coordinates=',xc, yc,slon, slat, xg, yg

    return
end subroutine Geo2Cart
