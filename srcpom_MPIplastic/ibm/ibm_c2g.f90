# 1 "ibm_c2g.F90"
subroutine Cart2Geo (xg, yg, slon, slat, xc, yc, pi, dg)
!! Converts Cartesian coordinates to Geographic coordinates

    real :: yg, xg          !Output fields: Geographic coordinates
    real :: slon, slat      !Starting Lon/Lat grid coordinates
    real :: xc, yc          !Cartesian coordinates
    real :: aa, pi, dg      !work field aa, pi=3.14..., dg=1degree in meters

    yg = slat + (yc / dg)
    aa = pi * yg / 180.
    xg = slon + (xc / (dg * cos(aa)))
!      print *, 'Geographic coordinates=', xg, yg

    return
end subroutine Cart2Geo
