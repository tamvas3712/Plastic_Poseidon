subroutine SIGMA2CART( Car, Sig, aland )
!! Convert Sigma coordinates to Cartesian.
    include 'comblk.h'
!    include 'grid.h'
    include 'fsm3d.h'
    real :: Car(im,jm,klev), Sig(im,jm,kb)

!   Assign values to last Sigma layer
    do j = 1, jm
        do i = 1, im
            Sig(i,j,kb) = Sig(i,j,kbm1)
        end do
    end do
!    write(*,*)'TESTSIG',kbm1,zz(1),h(15, 30)

!   for each grid point
    do 10 j = 1, jm
        do 15 i = 1, im
        !  check for land
            if(h(i, j) <= 1.) then
                do kk = 1, klev
                    Car(i, j, kk) = aland
                end do
                go to 15
            end if

        !   for each of std depth
            do 20 kc = 1, klev
            !  check if std depth greater than actual depth
                if(dpt(kc) > h(i,j)) then
                    Car(i, j, kc) = aland
                    go to 20
                end if
            !  find normalized depth
                sdpt = (-dpt(kc)) / h(i, j)
            !    if(i.eq.15.and.j.eq.30)write(*,*)'TESTSIG',kbm1,sdpt,zz(k)
            !  check if above first model's layer
                if(sdpt > zz(1)) then
                    Car(i, j, kc) = Sig(i, j, 1)
                    go to 20
                end if
            !  find appropriate layers for the std depth
            !  linear interpolation
                do 30 k = 1, kbm1
                    if(sdpt <= zz(k) .AND. sdpt >= zz(k+1)) then
                        Car(i, j, kc) = Sig(i, j, k+1) + &
                        (Sig(i, j, k) - Sig(i, j, k+1)) * &
                        (sdpt - zz(k+1)) / (zz(k) - zz(k+1))
                        go to 20
                    end if
!               if(i.eq.15.and.j.eq.30)write(*,*)'TESTSIG',kbm1 &
!                    kc,-dpt(kc),-zz(k)*h(i,j),Car(i,j,k),Sig(i,j,k), Sig(i, j, k+1)

                30 enddo
            20 enddo
        15 enddo
    10 enddo

    return
end subroutine SIGMA2CART
