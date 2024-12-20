# 1 "ibm_fcoord.F90"
subroutine FCoord (ip,jp,slon,slat,stepx,stepy,xc,yc,pi,dg)
!! Transform from Cartesian to Geographic and from Geographic to Grid coordinates
    real :: slon, slat, step, xc, yk, xg, yg

!      print *,'FCoord(Cart): xc,yc=',xc,yc

! Transform from Cartesian to Geographic coordinates
    call Cart2Geo (xg, yg, slon, slat, xc, yc, pi, dg)
!      print *,'Fcoord(Geo): xg,yg=',xg,yg

!     Transform from Geographic to Grid coordinates
    call Geo2Grid (ip,jp,slon,slat,stepx,stepy,xg,yg)
!      print *,'Fcoord(Grid): ip,jp=',ip,jp

    return
end subroutine FCoord
