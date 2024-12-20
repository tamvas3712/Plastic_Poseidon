# 1 "ibm_vertval.F90"
subroutine vertical_val(floc,i,j,im,jm,km,zz,f_sig,h,zp)
!! This routine returns the local vertical value on a given depth of a 3dfunction using linear interpolation.
!! 3d arrays are provided from POM and are in sigma coord.

    real :: f_sig(im, jm, km), h(im, jm)
    real :: zz(km)

    if(abs(zp) > h(i, j)) then
        floc = f_sig(i, j, km)
        return
    end if
!    print *, 'zp',zp
!    print *, zz

!   if particle is in the first half box
    if(zp > zz(1) * h(i, j)) then
        floc = f_sig(i, j, 1)
!        print *, 'first half ', floc
        go to 10
    end if
    do k = 2, km
        if(zp <= zz(k-1) * h(i, j) .AND. &
        zp > zz(k) * h(i, j)) then
            f1 = f_sig(i, j, k-1)
            f2 = f_sig(i, j, k)
            floc = f2 + (f1 - f2) * (zp - zz(k) * h(i, j)) / &
            (zz(k-1) * h(i, j) - zz(k) * h(i, j))
!            print *, i, j, k, f_sig(i, j, k), f_sig(i, j, k-1)
!            print *, 'floc',k,f1,f2,floc,zz(k-1)*h(i,j),zz(k)*h(i,j)
            go to 10
        end if
    end do
    10 continue

    return
end subroutine vertical_val
