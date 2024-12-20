# 1 "ibm_beach.F90"
!      subroutine beaching(rbt, ref, in, jn, bc, lbc, im, jm)
subroutine beaching(rbt, ref, in, jn, rtime, lbc, im, jm)
!!      This routine provides the retention time to the beach
!!      and a reflection flag which defines the way that the
!!      particle has to go back to the sea.
    integer :: ref, lbc(im, jm)
!      real    bc(im, jm), rbt
    real :: rtime,rbt
    double precision :: ro

!     Find retention time to the beach. Assume that the particle remaining at
!     the beach follows an exponential distribution.
    ro = rand(0)
!    ok = -alog(0.5) / bc(in, jn)
    ok = -alog(0.5) / rtime

    do 10 j = 1, 100
        if(ro > exp(-ok * j)) then
            rbt = float(j)
            go to 20
        end if
    10 END DO
    20 continue

! Assign a reflection flag according to the way the particle hits the beach
    ref = lbc(in, jn)
    if(lbc(in, jn) <= 0) then
        print *, 'error in ref value occured (not beached)'
    endif

    return
end subroutine beaching

