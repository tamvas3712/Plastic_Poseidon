subroutine INTERP_SIGMA(uloc,vloc,i,j,im,jm,km,dx,dy,zz,u,v,h,xp,yp,zp)
!! This routine provides local horizontal velocities from 3d arrays by linear interpolation.
!! Arrays are provided from POM and are in sigma coord.

    real :: u(im, jm, km), v(im, jm, km), h(im, jm)
    real :: zz(km)

    if(abs(zp) > h(i, j)) then
        uloc = 0.
        vloc = 0.
        return
    end if

!    print *, 'zp',zp
!    print *, zz

!   find u at left side
    hm = 0.5 * (h(i, j) + h(i-1, j))

!   if particle is in the first half box
    if(zp > zz(1) * hm) then
        ul = u(i, j, 1)
    !         print *, 'first half ', ul
        go to 10
    end if
    do k = 2, km
        if(zp <= zz(k-1) * hm .AND. &
        zp > zz(k) * hm) then
            u1 = u(i, j, k-1)
            u2 = u(i, j, k)
            ul = u2 + (u1 - u2) * (zp - zz(k) * hm) / &
            (zz(k-1) * hm - zz(k) * hm)
!            print *, i, j, k, u(i, j, k), u(i, j, k-1)
!            print *, 'ul',  k, u1, u2, ul, zz(k-1)*hm, zz(k)*hm
            go to 10
        end if
    end do
    10 continue

!   find u at right side
    hm = 0.5 * (h(i, j) + h(i+1, j))

!   if particle is in the first half box
    if(zp > zz(1) * hm) then
        ur = u(i+1, j, 1)
!        print *, 'first half ', ur
        go to 20
    end if
    do k = 2, km
        if(zp <= zz(k-1) * hm .AND. &
        zp > zz(k) * hm) then
            u1 = u(i+1, j, k-1)
            u2 = u(i+1, j, k)
            ur = u2 + (u1 - u2) * (zp - zz(k) * hm) / &
            (zz(k-1) * hm - zz(k) * hm)
!            print*, 'sig',  k, u1, u2, ur, zz(k-1)*hm, zz(k)*hm
            go to 20
        end if
    end do
    20 continue

!   find v at down side
    hm = 0.5 * (h(i, j) + h(i, j-1))

!   if particle is in the first half box
    if(zp > zz(1) * hm) then
        vd = v(i, j, 1)
!        print *, 'first half ', vd
        go to 30
    end if
    do k = 2, km
        if(zp <= zz(k-1) * hm .AND. zp > zz(k) * hm) then
            v1 = v(i, j, k-1)
            v2 = v(i, j, k)
            vd = v2 + (v1 - v2) * (zp - zz(k) * hm) / (zz(k-1) * hm - zz(k) * hm)
!            print *, i, j, k, v(i, j, k), v(i, j, k-1)
!            print *, 'vd',  k, v1, v2, vd, zz(k-1)*hm, zz(k)*hm
            go to 30
        end if
    end do
    30 continue

!   find v at upper side
    hm = 0.5 * (h(i, j) + h(i, j+1))

!   if particle is in the first half box
    if(zp > zz(1) * hm) then
        vu = v(i, j+1, 1)
!        print *, 'first half ', vd
        go to 40
    end if
    do k = 2, km
        if(zp <= zz(k-1) * hm .AND. &
        zp > zz(k) * hm) then
            v1 = v(i, j+1, k-1)
            v2 = v(i, j+1, k)
            vu = v2 + (v1 - v2) * (zp - zz(k) * hm) / &
            (zz(k-1) * hm - zz(k) * hm)
        !            print *, 'vu',  k, v1, v2, vu, zz(k-1)*hm, zz(k)*hm
            go to 40
        end if
    end do
    40 continue

    uloc = ul + (ur - ul) * (xp - (i - 1) * dx) / dx
    vloc = vd + (vu - vd) * (yp - (j - 1) * dy) / dy
!    print *, 'uloc, vloc, xp, yp, i, j', uloc, vloc, xp, yp, i, j


    return
end subroutine INTERP_SIGMA




subroutine INTERP_SIGMAav(uloc,vloc,i,j,im,jm,km,dx,dy,z,u,v,h,xp,yp,zp1,zp2)
!! This routine provides local horizontal velocities from 3d arrays by linear interpolation.
!! Arrays are provided from POM and are in sigma coord.

    real :: u(im, jm, km), v(im, jm, km), h(im, jm)
    real :: z(km),dz,ztot,zp1,zp2

! dim
!      if(abs(zp) .gt. h(i, j)) then
!         uloc = 0.
!         vloc = 0.
!         return
!      end if

!      print *, 'zp',zp
!      print *, z


!   find u at left and right side
    hm = 0.5*(h(i, j)+h(i+1,j))
    ztot=0.
    do k = 1, km
        if(ztot <= zp2 .AND. ztot >= zp1 )then
            if(k /= km)then
                dz=(z(k)-z(k+1))*hm
            else
                dz=(z(km-1)-z(km))*hm
            endif

            u1 = u(i, j, k)
            u2 = u(i+1, j, k)

            ul = ul+ u1 * dz
            ur = ur+ u2 * dz
            ztot=ztot+ dz
!            print *, i, j, k, u(i, j, k), u(i, j, k-1)
!            write(*,*) 'sigmaverage',  k, u1, u2, dz,ul
        end if
    end do
    ul = ul/ztot
    ur = ur/ztot
!      write(*,*) 'ul', ul,ur,ztot


!     find v at up and down side
    hm = 0.5 * (h(i, j) + h(i, j-1))
    ztot=0.
    do k = 1, km
        if(ztot <= zp2 .AND. ztot >= zp1 )then
            if(k /= km)then
                dz=(z(k)-z(k+1))*hm
            else
                dz=(z(km-1)-z(km))*hm
            endif
            v1 = v(i, j, k)
            v2 = v(i, j-1, k)
            vu = vu+ v1 * dz
            vd = vd+ v2 * dz
            ztot=ztot+ dz
        end if
    end do
    vd = vd/ztot
    vu = vu/ztot
    uloc = ul + (ur - ul) * (xp - (i - 1) * dx) / dx
    vloc = vd + (vu - vd) * (yp - (j - 1) * dy) / dy
!    print *, 'uloc', uloc, vloc, xp, yp, i, j

    return
end subroutine INTERP_SIGMAav

