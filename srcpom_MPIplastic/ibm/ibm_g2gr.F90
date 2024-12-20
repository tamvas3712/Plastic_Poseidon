subroutine Geo2Grid (ip,jp,slon,slat,stepx,stepy,xg,yg)
!! Transform from Geographic to Grid coordinates
    integer :: ip, jp       !Output fields: Particle Grid coordinates
    real ::    slon, slat   !Starting grid Geographic coordinates
    real ::    stepx, stepy !Model grid step
    real ::    xg, yg       !Particle Geographic coordinates

!     Transform from Geographic to Grid coordinates
    ip = int((xg-slon)/stepx)+1
    jp = int((yg-slat)/stepy)+1

!      print *,'Grid coordinates=', ip, jp

    return
end subroutine Geo2Grid
