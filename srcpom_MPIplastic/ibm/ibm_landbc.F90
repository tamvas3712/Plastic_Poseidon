!   land_bc.f    -- gt (c) 27/01/00

subroutine land_bc(lbc, h, im, jm, land_val)
!! this subroutine sets coastal boundary indices for particle reflections
!! It is assumed that=
!! Index -2 belongs to boxes outside of the field.      (escaped)
!! Index -1 belongs to bottom boxes.                    (bottom)
!! index 0 belongs to water boxes without land contact. (water)
!! index 1 belongs to land boxes without water contact. (land)
!! Index 2 to land boxes with water contact             (||)
!! to the left or to the  right.
!! Index 3 to land boxes with water contact             (=)
!! to the or to the bottom.
!! Index 4 to land boxes with more than                 (|_)
!! one direction of water box contact

    real ::    h(im, jm), land_val
    integer :: lbc(im, jm), iwork(im, jm)

    do j = 1, jm
        do i = 1,im
            iwork(i, j) = 0
            lbc(i, j) = 0
        end do
    end do

    do 10 j = 1, jm
        do 15 i = 1, im
            if(h(i, j) <= land_val) then
                iwork(i, j) = 1
            else
                iwork(i, j) = 0
            end if
        15 enddo
    10 enddo

    do 20 j = 2, jm - 1
        do 25 i = 2, im - 1
            ic = iwork(i, j)
            iu = iwork(i, j+1)
            id = iwork(i, j-1)
            il = iwork(i-1, j)
            ir = iwork(i+1, j)

            if(ic == 1) then
            !   for left and right
                if(ir == 0 ) then
                    lbc(i, j) = 2
                    goto 30
                end if
                if(il == 0) then
                    lbc(i, j) = 2
                    goto 30
                end if
                               
            !   for top and bottom
                if(iu == 0) then
                    lbc(i, j) = 3
                    goto 30
                end if
                if(id == 0) then
                    lbc(i, j) = 3
                    goto 30
                end if
                30 continue

        !     for land boxes adjacent with two water boxes making a corner
                if(ir == 0 .AND. iu == 0) then
                    lbc(i, j) = 4
                    goto 25
                end if
                if(ir == 0 .AND. id == 0) then
                    lbc(i, j) = 4
                    goto 25
                end if
                if(il == 0 .AND. iu == 0) then
                    lbc(i, j) = 4
                    goto 25
                end if
                if(il == 0 .AND. id == 0) then
                    lbc(i, j) = 4
                    goto 25
                end if
            end if
        25 enddo
    20 enddo

!   consider the limits
    do j = 1, jm
    !   left
        if(iwork(1, j) == 1) then
            if(iwork(2, j) == 0) lbc(1, j) = 2
            if(iwork(2, j) == 0 .AND. iwork(1, j+1) == 0) &
            lbc(1, j) = 4
            if(iwork(2, j) == 0 .AND. iwork(1, j-1) == 0) &
            lbc(1, j) = 4
        end if
    !   right
        if(iwork(im, j) == 1) then
            if(iwork(im-1, j) == 0) lbc(im, j) = 2
            if(iwork(im-1, j) == 0 .AND. iwork(im, j+1) == 0) &
            lbc(im, j) = 4
            if(iwork(im-1, j) == 0 .AND. iwork(im, j-1) == 0) &
            lbc(im, j) = 4
        end if
    end do

    do i = 1, im
    !   bottom
        if(iwork(i, 1) == 1) then
            if(iwork(i, 2) == 0) lbc(i, 1) = 3
            if(iwork(i, 2) == 0 .AND. iwork(i+1, 1) == 0) &
            lbc(i, 1) = 4
            if(iwork(i, 2) == 0 .AND. iwork(i-1, 1) == 0) &
            lbc(i, 1) = 4
        end if
    !   top
        if(iwork(i, jm) == 1) then
            if(iwork(i, jm-1) == 0) lbc(i, jm) = 3
            if(iwork(i, jm-1) == 0 .AND. iwork(i+1, jm) == 0) &
            lbc(i, jm) = 4
            if(iwork(i, jm-1) == 0 .AND. iwork(i-1, jm) == 0) &
            lbc(i, jm) = 4
        end if
    end do

    do j = 1, jm
        do i = 1,im
            if(lbc(i, j) == 0) lbc(i, j) = iwork(i, j)
        end do
    end do

    return
end subroutine land_bc
