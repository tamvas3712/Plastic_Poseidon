subroutine  goback(rbt,ref,xn,yn,xo,yo,dt,io,jo,im,jm,km,dx,dy,zz_os,u_os,v_os,h)
!! In case a particle is on the beach, find remaining beaching time (rbt).
!! If it has to go back to the sea, assign new coordinates taking into
!! account the way of beaching (ref).
!! Index -2 belongs to boxes outside of the field. (escaped)
!! Index -1 belongs to bottom boxes. (bottom)
!! Index 0 belongs to water boxes without land contact. (water)
!! Index 1 belongs to land boxes without water contact. (land)
!! Index 2 to land boxes with water contact to the left or to the  right.(||)
!! Index 3 to land boxes with water contact to the top or to the bottom.(=)
!! Index 4 to land boxes with more than one orientation to water box contact.(|_) 

    integer :: ref
    real :: u_os(im, jm, km), v_os(im, jm, km), h(im, jm)
    real :: zz_os(km)

    rbt = rbt - dt / 3600
!      write(*,*)'TEST REF',ref,rbt

    if(rbt <= 0) then
        call  INTERP_SIGMA(uloc,vloc,io,jo,im,jm,km,dx,dy,zz_os,u_os,v_os,h,xo,yo,0.)
!        print *, 'IO, JO', io, jo, ' ** REF', ref
        if(ref == 1) then
            xn = xo
            yn = yo
            rbt = 0.
            ref = 0
        else if(ref == 2) then
!            write(*,*)'GOBACK-2',xn,xo,uloc*dt,vloc * dt
            xn = xo
            yn = yo + vloc * dt
            rbt = 0.
            ref = 0
        else if(ref == 3) then
!            write(*,*)'GOBACK-3',xn,xo,uloc*dt,vloc * dt
            xn = xo + uloc * dt
            yn = yo
            rbt = 0.
            ref = 0
!          print *, 'IO, JO', io, jo, yn,yo,xn,xo
        else if(ref == 4) then
!            write(*,*)'GOBACK-4',xn,xo,uloc*dt,vloc * dt
            difx = abs(xn -xo)
            dify = abs(yn -yo)
            if(difx >= dify) then
                xn = xo + uloc * dt
                yn = yo
                rbt = 0.
                ref = 0
            else
                xn = xo
                yn = yo + vloc * dt
                rbt = 0.
                ref = 0
            endif
        else if(ref == 77) then
            xn = xo
            yn = yo
            rbt = 0.
            ref = 0
        else
            print *, 'error in particle beaching -  bad flag ref'
        end if
    end if

    return
end subroutine 
