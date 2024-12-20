! bounds_forcing.f

! spcify variable boundary conditions, atmospheric forcing, restoring

! ______________________________________________________________________
    subroutine bcond(idx)
!!  Apply open boundary conditions.
!!  Closed boundary conditions are automatically enabled through
!!  specification of the masks, dum, dvm and fsm, in which case the open
!!  boundary conditions, included below, will be overwritten
    implicit none
    include 'pom.h'
    integer :: idx
    integer :: i,j,k
    real :: ga,u1,wm
    real :: hmax
    real :: fxa,fxb,umid,vmid,denom,cl,vfac,ufac

    if(idx == 1) then

    ! eExternal (2-D) elevation boundary conditions
        do j=1,jm
            if(n_west == -1) elf(1,j)=elf(2,j)
            if(n_east == -1) elf(im,j)=elf(imm1,j)
        end do

        do i=1,im
            if(n_south == -1) elf(i,1)=elf(i,2)
            if(n_north == -1) elf(i,jm)=elf(i,jmm1)
        end do

        do j=1,jm
            do i=1,im
                elf(i,j)=elf(i,j)*fsm(i,j)
            end do
        end do

        return

    else if(idx == 2) then

    ! external (2-D) velocity boundary conditions
        do j=2,jmm1
        ! west
            if(n_west == -1) then
                uaf(2,j)=uabw(j)-rfw*sqrt(grav/h(2,j))*(el(2,j)-elw(j))
                uaf(2,j)=ramp*uaf(2,j)
                uaf(1,j)=uaf(2,j)
            !            vaf(1,j)=0.e0       !KTSIARAS
                umid=.5E0*(ua(2,j)+ua(2,j-1))
                vaf(1,j)=va(1,J)-dte/(dx(1,j)+dx(2,j)) &
                *(  (umid+abs(umid))*(va(1,j)-vabw(j)) &
                +(umid-ABS(umid))*(va(2,j)-va(1,j)) )

            end if

        ! east
            if(n_east == -1) then
                uaf(im,j)=uabe(j) &
                +rfe*sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
                uaf(im,j)=ramp*uaf(im,j)
            !           vaf(im,j)=0.e0        !KTSIARAS
                umid=.5E0*(ua(im,j)+ua(im,j-1))
                vaf(im,j)=va(im,j)-dte/(dx(im,j)+dx(imm1,j)) &
                *(  (umid+abs(umid))*(va(im,j)-va(imm1,j)) &
                +(umid-abs(umid))*(vabe(j) -va(im,j))  )

            end if
        end do

        do i=2,imm1
        ! south
            if(n_south == -1) then
                vaf(i,2)=vabs(i) &
                -rfs*sqrt(grav/h(i,2))*(el(i,2)-els(i))
                vaf(i,2)=ramp*vaf(i,2)
                vaf(i,1)=vaf(i,2)
            !            uaf(i,1)=0.e0        !KTSIARAS
                vmid=.5E0*(va(i,2)+va(i-1,2))
                uaf(i,1)=ua(i,1)-dte/(dy(i,1)+dy(i,2)) &
                *(   (vmid+abs(vmid))*(ua(i,1)-uabs(i)) &
                +(vmid-abs(vmid))*(ua(i,2)-ua(i,1)) )

            end if

        ! north
            if(n_north == -1) then
                vaf(i,jm)=vabn(i) &
                +rfn*sqrt(grav/h(i,jmm1))*(el(i,jmm1)-eln(i))
                vaf(i,jm)=ramp*vaf(i,jm)
            !           uaf(i,jm)=0.e0        !KTSIARAS
                vmid=.5E0*(va(i,jm)+va(i-1,jm))
                uaf(i,jm)=ua(i,jm)-dte/(dy(i,jm)+dy(i,jmm1)) &
                *(   (vmid+abs(vmid))*(ua(i,jm)-ua(i,jmm1)) &
                +(vmid-abs(vmid))*(uabn(i)-ua(i,jm))  )

            end if
        end do

        do j=1,jm
            do i=1,im
                uaf(i,j)=uaf(i,j)*dum(i,j)
                vaf(i,j)=vaf(i,j)*dvm(i,j)
            end do
        end do

    !        write(*,*)'TEST-BCOND-UVA',uaf(30,1),vaf(30,1),els(30),uabs(30)

        return

    else if(idx == 3) then

    ! internal (3-D) velocity boundary conditions
    ! radiation conditions
    ! smoothing is used in the direction tangential to the boundaries
        #ifdef orlanski
        hmax=4500.e0

        do k=1,kbm1
            do j=2,jmm1
            ! east
                if(n_east == -1) then
                    ga=sqrt(h(im,j)/hmax)
                    uf(im,j,k)=ga*(.25e0*u(imm1,j-1,k)+.5e0*u(imm1,j,k) &
                    +.25e0*u(imm1,j+1,k)) &
                    +(1.e0-ga)*(.25e0*u(im,j-1,k)+.5e0*u(im,j,k) &
                    +.25e0*u(im,j+1,k))
                    vf(im,j,k)=0.e0
                end if
            ! west
                if(n_west == -1) then
                    ga=sqrt(h(1,j)/hmax)
                    uf(2,j,k)=ga*(.25e0*u(3,j-1,k)+.5e0*u(3,j,k) &
                    +.25e0*u(3,j+1,k)) &
                    +(1.e0-ga)*(.25e0*u(2,j-1,k)+.5e0*u(2,j,k) &
                    +.25e0*u(2,j+1,k))
                    uf(1,j,k)=uf(2,j,k)
                    vf(1,j,k)=0.e0
                end if
            end do
        end do

        do k=1,kbm1
            do i=2,imm1
            ! south
                if(n_south == -1) then
                    ga=sqrt(h(i,1)/hmax)
                    vf(i,2,k)=ga*(.25e0*v(i-1,3,k)+.5e0*v(i,3,k) &
                    +.25e0*v(i+1,3,k)) &
                    +(1.e0-ga)*(.25e0*v(i-1,2,k)+.5e0*v(i,2,k) &
                    +.25e0*v(i+1,2,k))
                    vf(i,1,k)=vf(i,2,k)
                    uf(i,1,k)=0.e0
                end if
            ! north
                if(n_north == -1) then
                    ga=sqrt(h(i,jm)/hmax)
                    vf(i,jm,k)=ga*(.25e0*v(i-1,jmm1,k)+.5e0*v(i,jmm1,k) &
                    +.25e0*v(i+1,jmm1,k)) &
                    +(1.e0-ga)*(.25e0*v(i-1,jm,k)+.5e0*v(i,jm,k) &
                    +.25e0*v(i+1,jm,k))
                    uf(i,jm,k)=0.e0
                end if
            end do
        end do
        #else
    !--KTSIARAS
        #ifdef radiation
        do  k = 1, kbm1
            do  i=1,im
            ! north
                if(n_north == -1) then
                    denom=(vf(i,jm-1,k)+vb(i,jm-1,k)-2.*v(i,jm-2,k))
                    if(denom == 0.)denom=0.01
                    cl=(vb(i,jm-1,k)-vf(i,jm-1,k))/denom
                    if(cl > 1.)cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    vf(i,jm,k)=(vb(i,jm,k)*(1.-cl)+2.*cl*v(i,jm-1,k))/(1.+cl)
                    if (cl == 0.) vf(i,jm,k)=vbn(i,k)

                    vmid=.5E0*(v(i,jm,k)+v(i-1,jm,k))
                    uf(i,jm,k)=u(i,jm,k)-dti/(dy(i,jm)+dy(i,jmm1)) & !  **** NORTH
                    *(   (vmid+abs(vmid))*(u(i,jm,k)-u(i,jmm1,k)) &
                    +(vmid-abs(vmid))*(ubn(i,k)-u(i,jm,k))  )
                endif
            ! south
                if(n_south == -1) then
                    denom=(vf(i,3,k)+vb(i,3,k)-2.*v(i,4,k))
                    if(denom == 0.)denom=0.01
                    cl=(vb(i,3,k)-vf(i,3,k))/denom
                    if(cl > 1.)cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    vf(i,2,k)=(vb(i,2,k)*(1.-cl)+2.*cl*v(i,3,k))/(1.+cl)
                    if (cl == 0.) vf(i,2,k)=vbs(i,k)

                    vf(i,1,k)=vf(i,2,k)

                    vmid=.5E0*(v(i,2,k)+v(i-1,2,k))               ! **** SOUTH
                    uf(i,1,k)=u(i,1,k)-dti/(dy(i,1)+dy(i,2)) &
                    *(   (vmid+abs(vmid))*(u(i,1,k)-ubs(i,k)) &
                    +(vmid-abs(vmid))*(u(i,2,k)-u(i,1,k)) )
                endif
            enddo
        enddo

    !        write(*,*)'TEST-BCOND-UV',uf(30,1,1),vf(30,1,1),
    !     & ubs(30,1),vbs(30,1)

        do k = 1, kbm1
            do j = 1, jm
            ! east
                if(n_east == -1) then
                    uf(im,j,k)=ube(j,k)
                    denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.*u(im-2,j,k))
                    if(denom == 0.)denom=0.01
                    cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
                    if(cl > 1.)cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    uf(im,j,k)=(ub(im,j,k)*(1.-cl)+2.*cl*u(im-1,j,k))/(1.+cl)
                    if (cl == 0.) uf(im,j,k) = ube(j,k)

                    umid=.5E0*(u(im,j,k)+u(im,j-1,k))
                    vf(im,j,k)=v(im,j,k)-dti/(dx(im,j)+dx(imm1,j)) &
                    *(  (umid+abs(umid))*(v(im,j,k)-v(imm1,j,k)) &
                    +(umid-abs(umid))*(vbe(j,k)-v(im,j,k))  )
                endif
            ! west
                if(n_west == -1) then
                    denom=(uf(3,j,k)+ub(3,j,k)-2.*u(4,j,k))
                    if(denom == 0.)denom=0.01
                    cl=(ub(3,j,k)-uf(3,j,k))/denom
                    if(cl > 1.)cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    uf(2,j,k)=(ub(2,j,k)*(1.-cl)+2.*cl*u(3,j,k))/(1.+cl)
                    if (cl == 0.) uf(2,j,k) = ubw(j,k)
                    uf(1,j,k)=uf(2,j,k)

                    umid=.5E0*(u(2,j,k)+u(2,j-1,k))
                    vf(1,j,k)=v(1,j,k)-dti/(dx(1,j)+dx(2,j)) &
                    *(  (umid+abs(umid))*(v(1,j,k)-vbw(j,k)) &
                    +(umid-abs(umid))*(v(2,j,k)-v(1,j,k)) )

                endif
            enddo
        enddo

        #else

        do  k = 1, kbm1
            do  i=1,im
            ! north
                if(n_north == -1) then
                    vf(i,jm,k)=vbn(i,k)
                    uf(i,jm,k)=ubn(i,k)
                endif
            ! south
                if(n_south == -1) then
                    vf(i,2,k)=vbs(i,k)
                    vf(i,1,k)=vf(i,2,k)
                    uf(i,1,k)=ubs(i,k)
                endif
            enddo
        enddo

        do k = 1, kbm1
            DO j = 1, jm
            ! east
                if(n_east == -1) then
                    uf(im,j,k)=ube(j,k)
                    vf(im,j,k)=vbe(j,k)
                
                endif
            ! west
                if(n_west == -1) then
                    uf(2,j,k)=ubw(j,k)
                    uf(1,j,k)=uf(2,j,k)
                    vf(1,j,k)=vbw(j,k)
                endif
            enddo
        enddo
        #endif
        #endif

        do k=1,kbm1
            do j=1,jm
                do i=1,im
                    uf(i,j,k)=uf(i,j,k)*dum(i,j)
                    vf(i,j,k)=vf(i,j,k)*dvm(i,j)
                end do
            end do
        end do

        return

    else if(idx == 4) then

    ! temperature and salinity boundary conditions (using uf and vf,
    ! respectively)

        #ifdef palma
        do k=1,kbm1
            do i=1,im
            ! north
                if(n_north == -1) then
                ! temperature
                    fxb=(uf(i,jm-1,k)+tb(i,jm-1,k)-2.*t(i,jm-2,k))
                    fxa=(tb(i,jm-1,k)-uf(i,jm-1,k))
                    vfac=2.*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
                    if (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        vfac=2.*(v(i,jm,k)+cl)*dti/(dy(i,jm)+dy(i,jmm1))
                    end if
                    if(vfac >= 0.) then
                        uf(i,jm,k)=t(i,jm,k)-vfac*(t(i,jm,k)-t(i,jmm1,k))
                    else
                        uf(i,jm,k)=t(i,jm,k)-vfac*(tbn(i,k)-t(i,jm,k))
                    endif
                ! salinity
                    fxb=(vf(i,jm-1,k)+sb(i,jm-1,k)-2.*s(i,jm-2,k))
                    fxa=(sb(i,jm-1,k)-vf(i,jm-1,k))
                    vfac=2.*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
                    if (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        vfac=2.*(v(i,jm,k)+cl)*dti/(dy(i,jm)+dy(i,jmm1))
                    end if
                    if(vfac >= 0.) then
                        vf(i,jm,k)=s(i,jm,k)-vfac*(s(i,jm,k)-s(i,jmm1,k))
                    else
                        vf(i,jm,k)=s(i,jm,k)-vfac*(sbn(i,k)-s(i,jm,k))
                    endif
                endif
            ! south
                if(n_south == -1) then
                ! temperature
                    fxb=(uf(i,2,k)+tb(i,2,k)-2.*t(i,3,k))
                    fxa=(tb(i,2,k)-uf(i,2,k))
                    vfac=2.*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
                    if (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        vfac=2.*(v(i,2,k)-cl)*dti/(dy(i,1)+dy(i,2))
                    end if
                    if(vfac <= 0.) then
                        uf(i,1,k)=t(i,1,k)-vfac*(t(i,2,k)-t(i,1,k))
                    else
                        uf(i,1,k)=t(i,1,k)-vfac*(t(i,1,k)-tbs(i,k))
                    endif
                ! salinity
                    fxb=(vf(i,2,k)+sb(i,2,k)-2.*s(i,3,k))
                    fxa=(sb(i,2,k)-vf(i,2,k))
                    vfac=2.*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
                    if (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        vfac=2.*(v(i,2,k)-cl)*dti/(dy(i,1)+dy(i,2))
                    end if
                    if(vfac <= 0.) then
                        vf(i,1,k)=s(i,1,k)-vfac*(s(i,2,k)-s(i,1,k))
                    else
                        vf(i,1,k)=s(i,1,k)-vfac*(s(i,1,k)-sbs(i,k))
                    endif
                endif
            enddo
        enddo

    !        write(*,*)'TEST-BCOND-TS',uf(30,1,1),vf(30,1,1),
    !     & tbs(30,1),sbs(30,1)

        do k=1,kbm1
            do j=1,jm
            ! east
                if(n_east == -1) then
                ! temperature
                    fxb=(uf(im-1,j,k)+tb(im-1,j,k)-2.*t(im-2,j,k))
                    fxa=(tb(im-1,j,k)-uf(im-1,j,k))
                    ufac=2.*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
                    if (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        ufac=2.*(u(im,j,k)+cl)*dti/(dx(im,j)+dx(imm1,j))
                    end if
                    if(ufac >= 0.) then
                        uf(im,j,k)=t(im,j,k)-ufac*(t(im,j,k)-t(imm1,j,k))
                    else
                        uf(im,j,k)=t(im,j,k)-ufac*(tbe(j,k)-t(im,j,k))
                    endif
                ! salinity
                    fxb=(vf(im-1,j,k)+sb(im-1,j,k)-2.*s(im-2,j,k))
                    fxa=(sb(im-1,j,k)-vf(im-1,j,k))
                    ufac=2.*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
                    if (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        ufac=2.*(u(im,j,k)+cl)*dti/(dx(im,j)+dx(imm1,j))
                    end if
                    if(ufac >= 0.) then
                        vf(im,j,k)=s(im,j,k)-ufac*(s(im,j,k)-s(imm1,j,k))
                    else
                        vf(im,j,k)=s(im,j,k)-ufac*(sbe(j,k)-s(im,j,k))
                    endif
                endif
            ! west
                if(n_west == -1) then
                ! temperature
                    fxb=(uf(2,j,k)+tb(2,j,k)-2.*t(3,j,k))
                    fxa=(tb(2,j,k)-uf(2,j,k))
                    ufac=2.*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
                    if (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        ufac=2.*(u(2,j,k)-cl)*dti/(dx(1,j)+dx(2,j))
                    end if
                    if(ufac <= 0.) then
                        uf(1,j,k)=t(1,j,k)-ufac*(t(2,j,k)-t(1,j,k))
                    else
                        uf(1,j,k)=t(1,j,k)-ufac*(t(1,j,k)-tbw(j,k))
                    endif
                ! salinity
                    fxb=(vf(2,j,k)+sb(2,j,k)-2.*s(3,j,k))
                    fxa=(sb(2,j,k)-vf(2,j,k))
                    ufac=2.*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
                    IF (fxa*fxb > 0. .AND. abs(fxb) > abs(fxa)) then
                        cl = fxa/fxb
                        ufac=2.*(u(2,j,k)-cl)*dti/(dx(1,j)+dx(2,j))
                    end if
                    if(ufac <= 0.) then
                        vf(1,j,k)=s(1,j,k)-ufac*(s(2,j,k)-s(1,j,k))
                    else
                        vf(1,j,k)=s(1,j,k)-ufac*(s(1,j,k)-sbw(j,k))
                    endif
                endif
            enddo
        enddo
        #else

        do k=1,kbm1
            do j=1,jm
            ! east
                if(n_east == -1) then
                    u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
                    if(u1 <= 0.e0) then
                        uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
                        vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
                    else
                        uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
                        vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
                        if(k /= 1 .AND. k /= kbm1) then
                            wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti &
                            /((zz(k-1)-zz(k+1))*dt(imm1,j))
                            uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                            vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
                        endif
                    end if
                end if

            ! west
                if(n_west == -1) then
                    u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
                    if(u1 >= 0.e0) then
                        uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
                        vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
                    else
                        uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
                        vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
                        if(k /= 1 .AND. k /= kbm1) then
                            wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti &
                            /((zz(k-1)-zz(k+1))*dt(2,j))
                            uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                            vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
                        end if
                    end if
                end if
            end do

            do i=1,im
            ! south
                if(n_south == -1) then
                    u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
                    if(u1 >= 0.e0) then
                        uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
                        vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
                    else
                        uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
                        vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
                        if(k /= 1 .AND. k /= kbm1) then
                            wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti &
                            /((zz(k-1)-zz(k+1))*dt(i,2))
                            uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                            vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
                        end if
                    end if
                end if

            ! north
                if(n_north == -1) then
                    u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
                    if(u1 <= 0.e0) then
                        uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
                        vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
                    else
                        uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
                        vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
                        if(k /= 1 .AND. k /= kbm1) then
                            wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti &
                            /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                            uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                            vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
                        end if
                    end if
                end if
            end do
        end do

        do k=1,kbm1
            do j=1,jm
                do i=1,im
                    uf(i,j,k)=uf(i,j,k)*fsm(i,j)
                    vf(i,j,k)=vf(i,j,k)*fsm(i,j)
                end do
            end do
        end do

        return

    else if(idx == 5) then

    ! vertical velocity boundary conditions
        do k=1,kbm1
            do j=1,jm
                do i=1,im
                    w(i,j,k)=w(i,j,k)*fsm(i,j)
                end do
            end do
        end do

        return

    else if(idx == 6) then

    ! q2 and q2l boundary conditions

        do k=1,kb
            do j=1,jm
            ! west
                if(n_west == -1) then
                    u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
                    if(u1 >= 0.e0) then
                        uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
                        vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
                    else
                        uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
                        vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
                    end if
                end if

            ! east
                if(n_east == -1) then
                    u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
                    if(u1 <= 0.e0) then
                        uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
                        vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
                    else
                        uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
                        vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
                    end if
                end if
            end do

            do i=1,im
            ! south
                if(n_south == -1) then
                    u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
                    if(u1 >= 0.e0) then
                        uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
                        vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
                    else
                        uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
                        vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
                    end if
                end if

            ! north
                if(n_north == -1) then
                    u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
                    if(u1 <= 0.e0) then
                        uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
                        vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
                    else
                        uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
                        vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
                    end if
                end if
            end do
        end do

        do k=1,kb
            do j=1,jm
                do i=1,im
                    uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
                    vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
                end do
            end do
        end do

        return

    endif

    end subroutine bcond

! ______________________________________________________________________
    subroutine bcondorl(idx)
!!  This is an optional subroutine replacing bcond (applying boundary conditions) and using Orlanski's
!!  scheme (J. Comp. Phys. 21, 251-269, 1976), specialized for the
!!  seamount problem
    implicit none
    include 'pom.h'
    integer :: idx
    real :: cl,denom
    real :: ar,eps
    integer :: i,j,k

    if(idx == 1) then

    ! external (2-D) elevation boundary conditions
        do  j=1,jm
            if(n_west == -1) elf(1,j)=elf(2,j)
            if(n_east == -1) elf(im,j)=elf(imm1,j)
        end do

        do j=1,jm
            do i=1,im
                elf(i,j)=elf(i,j)*fsm(i,j)
            end do
        end do

        return

    else if(idx == 2) then

    ! external (2-D) velocity  boundary conditions
        do j=2,jmm1
        ! east
            if(n_east == -1) then
                denom=(uaf(im-1,j)+uab(im-1,j)-2.e0*ua(im-2,j))
                if(denom == 0.0e0)denom=0.01e0
                cl=(uab(im-1,j)-uaf(im-1,j))/denom
                if(cl > 1.e0) cl=1.e0
                if(cl < 0.e0) cl=0.e0
                uaf(im,j)=(uab(im,j)*(1.e0-cl)+2.e0*cl*ua(im-1,j)) &
                /(1.e0+cl)
                vaf(im,j)=0.e0
            end if

        ! west
            if(n_west == -1) then
                denom=(uaf(3,j)+uab(3,j)-2.e0*ua(4,j))
                if(denom == 0.0e0)denom=0.01e0
                cl=(uab(3,j)-uaf(3,j))/denom
                if(cl > 1.e0) cl=1.e0
                if(cl < 0.e0) cl=0.e0
                uaf(2,j)=(uab(2,j)*(1.e0-cl)+2.e0*cl*ua(3,j)) &
                /(1.e0+cl)
                uaf(1,j)=uaf(2,j)
                vaf(1,j)=0.e0
            end if
        end do

        do i=2,imm1
        ! south
            if(n_south == -1) then
                denom=(vaf(i,3)+vab(i,3)-2.e0*va(i,4))
                if(denom == 0.0e0)denom=0.01e0
                cl=(vab(i,3)-vaf(i,3))/denom
                if(cl > 1.e0) cl=1.e0
                if(cl < 0.e0) cl=0.e0
                vaf(i,2)=(vab(i,2)*(1.e0-cl)+2.e0*cl*va(i,3)) &
                /(1.e0+cl)
                vaf(i,1)=vaf(i,2)
                uaf(i,1)=0.e0
            end if

        ! north
            if(n_north == -1) then
                denom=(vaf(i,jm-1)+vab(i,jm-1)-2.e0*va(i,jm-2))
                if(denom == 0.0e0)denom=0.01e0
                cl=(vab(i,jm-1)-vaf(i,jm-1))/denom
                if(cl > 1.e0) cl=1.e0
                if(cl < 0.e0) cl=0.e0
                vaf(i,jm)=(vab(i,jm)*(1.e0-cl)+2.e0*cl*va(i,jm-1)) &
                /(1.e0+cl)
                uaf(i,jm)=0.e0
            end if
        end do

        do j=1,jm
            do i=1,im
                uaf(i,j)=uaf(i,j)*dum(i,j)
                vaf(i,j)=vaf(i,j)*dvm(i,j)
            end do
        end do

        return

    else if(idx == 3) then

    ! internal (3-D) velocity boundary conditions

        do k=1,kbm1
            do j=2,jmm1
            ! east
                if(n_east == -1) then
                    denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.e0*u(im-2,j,k))
                    if(denom == 0.e0)denom=0.01e0
                    cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    uf(im,j,k)=(ub(im,j,k)*(1.e0-cl)+2.e0*cl*u(im-1,j,k)) &
                    /(1.e0+cl)
                    vf(im,j,k)=0.e0
                end if

            ! west
                if(n_west == -1) then
                    denom=(uf(3,j,k)+ub(3,j,k)-2.e0*u(4,j,k))
                    if(denom == 0.e0)denom=0.01e0
                    cl=(ub(3,j,k)-uf(3,j,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    uf(2,j,k)=(ub(2,j,k)*(1.e0-cl)+2.e0*cl*u(3,j,k)) &
                    /(1.e0+cl)
                    uf(1,j,k)=uf(2,j,k)
                    vf(1,j,k)=0.e0
                end if
            end do

            do i=2,imm1
            ! south
                if(n_south == -1) then
                    denom=(vf(i,3,k)+vb(i,3,k)-2.e0*v(i,4,k))
                    if(abs(denom) == 0.0e0)denom=0.01e0
                    cl=(vb(i,3,k)-vf(i,3,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    vf(i,2,k)=(vb(i,2,k)*(1.e0-cl)+2.e0*cl*v(i,3,k)) &
                    /(1.e0+cl)
                    vf(i,1,k)=vf(i,2,k)
                    uf(i,1,k)=0.e0
                end if

            ! north
                if(n_north == -1) then
                    denom=(vf(i,jm-1,k)+vb(i,jm-1,k)-2.e0*v(i,jm-2,k))
                    if(abs(denom) == 0.0e0)denom=0.01e0
                    cl=(vb(i,jm-1,k)-vf(i,jm-1,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    vf(i,jm,k)=(vb(i,jm,k)*(1.e0-cl)+2.e0*cl*v(i,jm-1,k)) &
                    /(1.e0+cl)
                    uf(i,jm,k)=0.e0
                end if
            end do
        end do

        do k=1,kbm1
            do j=1,jm
                do i=1,im
                    uf(i,j,k)=uf(i,j,k)*dum(i,j)
                    vf(i,j,k)=vf(i,j,k)*dvm(i,j)
                end do
            end do
        end do

        return

    else if(idx == 4) then

    ! temperature and salinity boundary conditions (using uf and vf,
    ! respectively)
        do k=1,kbm1
            do i=1,im
            ! north KTSIARAS
                if(n_north == -1) then
                    denom=(uf(i,jm-1,k)+tb(i,jm-1,k)-2.*t(i,jm-2,k))
                    if(denom == 0.)denom=0.01
                    cl=(tb(i,jm-1,k)-uf(i,jm-1,k))/denom
                    if(cl > 1.) cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    uf(i,jm,k)=(tb(i,jm,k)*(1.-cl)+2.*cl*t(i,jm-1,k))/(1.+cl)
                    if(cl == 0. .AND. vb(i,jm,k) <= 0.) uf(i,jm,k)=tbn(i,k)
                
                    denom=(vf(i,jm-1,k)+sb(i,jm-1,k)-2.*s(i,jm-2,k))
                    if(denom == 0.)denom=0.01
                    cl=(sb(i,jm-1,k)-vf(i,jm-1,k))/denom
                    if(cl > 1.) cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    vf(i,jm,k)=(sb(i,jm,k)*(1.-cl)+2.*cl*s(i,jm-1,k))/(1.+cl)
                    if(cl == 0. .AND. vb(i,jm,k) <= 0.) vf(i,jm,k)=sbn(i,k)
                endif
            ! south KTSIARAS
                if(n_south == -1) then
                    denom=(uf(i,2,k)+tb(i,2,k)-2.*t(i,3,k))
                    if(denom == 0.)denom=0.01
                    cl=(tb(i,2,k)-uf(i,2,k))/denom
                    if(cl > 1.) cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    uf(i,1,k)=(tb(i,1,k)*(1.-cl)+2.*cl*t(i,2,k))/(1.+cl)
                    if(cl == 0. .AND. vb(i,2,k) >= 0.) uf(i,1,k)=tbs(i,k)
                
                    denom=(vf(i,2,k)+sb(i,2,k)-2.*s(i,3,k))
                    if(denom == 0.)denom=0.01
                    cl=(sb(i,2,k)-vf(i,2,k))/denom
                    if(cl > 1.) cl=1.
                    if(cl < 0.) then
                        cl=0.
                    endif
                    vf(i,1,k)=(sb(i,1,k)*(1.-cl)+2.*cl*s(i,2,k))/(1.+cl)
                    if(cl == 0. .AND. vb(i,2,k) >= 0.) vf(i,1,k)=sbs(i,k)
                endif
            enddo
        enddo

        do k=1,kbm1
            do j=1,jm
            ! east
                if(n_east == -1) then
                    ube(j,k)=ub(im,j,k)
                    denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.e0*t(im-2,j,k))
                    if(denom == 0.e0) denom=0.01e0
                    cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    uf(im,j,k)=(tb(im,j,k)*(1.e0-cl)+2.e0*cl*t(im-1,j,k)) &
                    /(1.e0+cl)
                    if(cl == 0.e0 .AND. ube(j,k) <= 0.e0) uf(im,j,k)=tbe(j,k)

                    denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.e0*s(im-2,j,k))
                    if(denom == 0.e0) denom=0.01e0
                    cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    vf(im,j,k)=(sb(im,j,k)*(1.e0-cl)+2.e0*cl*s(im-1,j,k)) &
                    /(1.e0+cl)
                    if(cl == 0.e0 .AND. ube(j,k) <= 0.e0) vf(im,j,k)=sbe(j,k)
                end if

            ! west
                if(n_west == -1) then
                    ubw(j,k)=ub(2,j,k)
                    denom=(uf(2,j,k)+tb(2,j,k)-2.e0*t(3,j,k))
                    if(denom == 0.e0) denom=0.01e0
                    cl=(tb(2,j,k)-uf(2,j,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    uf(1,j,k)=(tb(1,j,k)*(1.e0-cl)+2.e0*cl*t(2,j,k))/(1.e0+cl)
                    if(cl == 0.e0 .AND. ubw(j,k) >= 0.e0) uf(1,j,k)=tbw(j,k)

                    denom=(vf(2,j,k)+sb(2,j,k)-2.e0*s(3,j,k))
                    if(denom == 0.e0) denom=0.01e0
                    cl=(sb(2,j,k)-vf(2,j,k))/denom
                    if(cl > 1.e0) cl=1.e0
                    if(cl < 0.e0) cl=0.e0
                    vf(1,j,k)=(sb(1,j,k)*(1.e0-cl)+2.e0*cl*s(2,j,k))/(1.e0+cl)
                    if(cl == 0.e0 .AND. ubw(j,k) >= 0.e0) vf(1,j,k)=sbw(j,k)
                end if
            end do
        end do

        do k=1,kbm1
            do j=1,jm
                do i=1,im
                    uf(i,j,k)=uf(i,j,k)*fsm(i,j)
                    vf(i,j,k)=vf(i,j,k)*fsm(i,j)
                end do
            end do
        end do

        return

    else if(idx == 5) then

    ! vertical velocity boundary conditions
        do k=1,kbm1
            do j=1,jm
                do i=1,im
                    w(i,j,k)=w(i,j,k)*fsm(i,j)
                end do
            end do
        end do

        do k=1,kb   !KTSIARAS
            do j=1,jm
                if(n_east == -1) then
                    w(im,j,k)=0.E0       ! **** EAST
                endif
                if(n_east == -1) then
                    w(1,j,k)=0.E0        ! **** WEST
                endif
            enddo
            do i=1,im
                if(n_north == -1) then
                    w(i,jm,k)=0.E0       ! **** NORTH
                endif
                if(n_south == -1) then
                    w(i,1,k)=0.E0        ! **** SOUTH
                endif
            enddo
        enddo

        return

    else if(idx == 6) then

    ! q2 and q2l boundary conditions
        do k=1,kb
            do j=1,jm
                if(n_east == -1) then
                    uf(im,j,k)=1.e-10
                    vf(im,j,k)=1.e-10
                end if
                if(n_west == -1) then
                    uf(1,j,k)=1.e-10
                    vf(1,j,k)=1.e-10
                end if
            end do

            do j=1,jm
                do i=1,im
                    uf(i,j,k)=uf(i,j,k)*fsm(i,j)
                    vf(i,j,k)=vf(i,j,k)*fsm(i,j)
                end do
            end do
        end do

        return

    endif

    end subroutine bcondorl

! _____________________________________________________________________
    subroutine lateral_bc
!!  Create variable lateral boundary conditions.
    implicit none
    include 'pom.h'
    integer :: nz
    parameter(nz=40)
    integer :: i,j,k,ntime,ibc
    real :: tbc,fold,fnew
    real :: z0(nz),hs(im,jm)
    real :: t_w(jm,nz),s_w(jm,nz),t_e(jm,nz),s_e(jm,nz),q, &
    t_n(im,nz),s_n(im,nz),t_s(im,nz),s_s(im,nz)
    real :: ts1(im,jm,nz),ss1(im,jm,nz)
    real :: ts2(im,jm,kb),ss2(im,jm,kb)

    tbc=30 ! time between bc files (days)
    ibc=int(tbc*86400.e0/dti)
    ntime=time/tbc
! read bc data
! read initial bc file
    if (iint == 1) then
        call read_boundary_conditions_pnetcdf(iint/ibc,nz,z0, &
        t_w,s_w,uabwf,t_e,s_e,uabef,t_n,s_n,vabnf,t_s,s_s,vabsf)
    !  south
        do i=1,im
            do j=1,jm
                hs(i,j)=h(i,2)
                do k=1,nz
                    ts1(i,j,k)=t_s(i,k)
                    ss1(i,j,k)=s_s(i,k)
                end do
            end do
        end do
        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        do i=1,im
            do k=1,kb
                tbsf(i,k)=ts2(i,jm/2,k)
                sbsf(i,k)=ss2(i,jm/2,k)
            end do
        end do
    ! north
        do i=1,im
            do j=1,jm
                hs(i,j)=h(i,jmm1)
                do k=1,nz
                    ts1(i,j,k)=t_n(i,k)
                    ss1(i,j,k)=s_n(i,k)
                end do
            end do
        end do
        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        do i=1,im
            do k=1,kb
                tbnf(i,k)=ts2(i,jm/2,k)
                sbnf(i,k)=ss2(i,jm/2,k)
            end do
        end do
    ! east
        do i=1,im
            do j=1,jm
                hs(i,j)=h(imm1,j)
                do k=1,nz
                    ts1(i,j,k)=t_e(j,k)
                    ss1(i,j,k)=s_e(j,k)
                end do
            end do
        end do
        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        do j=1,jm
            do k=1,kb
                tbef(j,k)=ts2(im/2,j,k)
                sbef(j,k)=ss2(im/2,j,k)
            end do
        end do
    ! west
        do i=1,im
            do j=1,jm
                hs(i,j)=h(2,j)
                do k=1,nz
                    ts1(i,j,k)=t_w(j,k)
                    ss1(i,j,k)=s_w(j,k)
                end do
            end do
        end do
        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        do j=1,jm
            do k=1,kb
                tbwf(j,k)=ts2(im/2,j,k)
                sbwf(j,k)=ss2(im/2,j,k)
            end do
        end do

    end if
! read bc file corresponding to next time
    if (iint == 1 .OR. mod(iint,ibc) == 0.) then
        do j=1,jm
            do k=1,kb
                tbwb(j,k)=tbwf(j,k)
                sbwb(j,k)=sbwf(j,k)
                tbeb(j,k)=tbef(j,k)
                sbeb(j,k)=sbef(j,k)
            end do
            uabwb(j)=uabwf(j)
            uabeb(j)=uabef(j)
        end do
        do i=1,im
            do k=1,kb
                tbnb(i,k)=tbnf(i,k)
                sbnb(i,k)=sbnf(i,k)
                tbsb(i,k)=tbsf(i,k)
                sbsb(i,k)=sbsf(i,k)
            end do
            vabnb(i)=vabnf(i)
            vabsb(i)=vabsf(i)
        end do
        if (iint /= iend) then
            call read_boundary_conditions_pnetcdf((iint+ibc)/ibc,nz,z0, &
            t_w,s_w,uabwf,t_e,s_e,uabef,t_n,s_n,vabnf,t_s,s_s,vabsf)
        ! south
            do i=1,im
                do j=1,jm
                    hs(i,j)=h(i,2)
                    do k=1,nz
                        ts1(i,j,k)=t_s(i,k)
                        ss1(i,j,k)=s_s(i,k)
                    end do
                end do
            end do
            call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            do i=1,im
                do k=1,kb
                    tbsf(i,k)=ts2(i,jm/2,k)
                    sbsf(i,k)=ss2(i,jm/2,k)
                end do
            end do
        ! north
            do i=1,im
                do j=1,jm
                    hs(i,j)=h(i,jmm1)
                    do k=1,nz
                        ts1(i,j,k)=t_n(i,k)
                        ss1(i,j,k)=s_n(i,k)
                    end do
                end do
            end do
            call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            do i=1,im
                do k=1,kb
                    tbnf(i,k)=ts2(i,jm/2,k)
                    sbnf(i,k)=ss2(i,jm/2,k)
                end do
            end do
        ! east
            do i=1,im
                do j=1,jm
                    hs(i,j)=h(imm1,j)
                    do k=1,nz
                        ts1(i,j,k)=t_e(j,k)
                        ss1(i,j,k)=s_e(j,k)
                    end do
                end do
            end do
            call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            do j=1,jm
                do k=1,kb
                    tbef(j,k)=ts2(im/2,j,k)
                    sbef(j,k)=ss2(im/2,j,k)
                end do
            end do
        ! west
            do i=1,im
                do j=1,jm
                    hs(i,j)=h(2,j)
                    do k=1,nz
                        ts1(i,j,k)=t_w(j,k)
                        ss1(i,j,k)=s_w(j,k)
                    end do
                end do
            end do
            call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            do j=1,jm
                do k=1,kb
                    tbwf(j,k)=ts2(im/2,j,k)
                    sbwf(j,k)=ss2(im/2,j,k)
                end do
            end do
        ! end
        end if
    end if

! linear interpolation in time
    fnew=time/tbc-ntime
    fold=1.-fnew
    do j=1,jm
        do k=1,kb
            tbw(j,k)=fold*tbwb(j,k)+fnew*tbwf(j,k)
            sbw(j,k)=fold*sbwb(j,k)+fnew*sbwf(j,k)
            tbe(j,k)=fold*tbeb(j,k)+fnew*tbef(j,k)
            sbe(j,k)=fold*sbeb(j,k)+fnew*sbef(j,k)
        end do
        uabe(j)=fold*uabeb(j)+fnew*uabef(j)
    end do
    do i=1,im
        do k=1,kb
            tbn(i,k)=fold*tbnb(i,k)+fnew*tbnf(i,k)
            sbn(i,k)=fold*sbnb(i,k)+fnew*sbnf(i,k)
            tbs(i,k)=fold*tbsb(i,k)+fnew*tbsf(i,k)
            sbs(i,k)=fold*sbsb(i,k)+fnew*sbsf(i,k)
        end do
        vabn(i)=fold*vabnb(i)+fnew*vabnf(i)
        vabs(i)=fold*vabsb(i)+fnew*vabsf(i)
    end do

    return
    end subroutine lateral_bc

! ______________________________________________________________________
    subroutine wind
!!  Read and interpolate (in time) wind stress.
    implicit none
    include 'pom.h'
    integer :: i,j,ntime,iwind
    real :: twind,fold,fnew

    twind=30 ! time between wind files (days)
    iwind=int(twind*86400.e0/dti)

! read wind stress data
! read initial wind file
    if (iint == 1) call read_wind_pnetcdf(iint/iwind,wusurff,wvsurff)
! read wind file corresponding to next time
    if (iint == 1 .OR. mod(iint,iwind) == 0.) then
        do i=1,im
            do j=1,jm
                wusurfb(i,j)=wusurff(i,j)
                wvsurfb(i,j)=wvsurff(i,j)
            end do
        end do
        if (iint /= iend) then
            call read_wind_pnetcdf((iint+iwind)/iwind,wusurff,wvsurff)
        end if
    end if

! linear interpolation in time
    ntime=time/twind
    fnew=time/twind-ntime
    fold=1.-fnew
    do i=1,im
        do j=1,jm
            wusurf(i,j)=fold*wusurfb(i,j)+fnew*wusurff(i,j)
            wvsurf(i,j)=fold*wvsurfb(i,j)+fnew*wvsurff(i,j)
        end do
    end do

    return
    end subroutine wind

! ______________________________________________________________________
    subroutine heat
!!  Read and interpolate heat flux in time
    implicit none
    include 'pom.h'
    integer :: i,j,ntime,iheat
    real :: theat,fold,fnew

    theat=1. ! time between wind forcing (days)
    iheat=int(theat*86400.e0/dti)

! read wind stress data
! read initial heat file
    if (iint == 1) call read_heat_pnetcdf(iint/iheat,wtsurff,swradf)
! read heat forcing corresponding to next twind
    if (iint == 1 .OR. mod(iint,iheat) == 0.) then
        do i=1,im
            do j=1,jm
                wtsurfb(i,j)=wtsurff(i,j)
                swradb(i,j)=swradf(i,j)
            end do
        end do
        call read_heat_pnetcdf((iint+iheat)/iheat,wtsurff,swradf)
    end if

! linear interpolation in time
    ntime=time/theat
    fnew=time/theat-ntime
    fold=1.-fnew
    do i=1,im
        do j=1,jm
            wtsurf(i,j)=fold*wtsurfb(i,j)+fnew*wtsurff(i,j)
            swrad(i,j)=fold*swradb(i,j)+fnew*swradf(i,j)
        end do
    end do

    return
    end subroutine heat

! ______________________________________________________________________
    subroutine water
!!  Read and interpolate water flux in time
    implicit none
    include 'pom.h'
    integer :: i,j,ntime,iwater
    real :: twater,fold,fnew

    twater=1. ! time between wind forcing (days)
    iwater=int(twater*86400.e0/dti)

! read wind stress data
! read initial water file
    if (iint == 1) call read_water_pnetcdf(iint/iwater,wssurff)
! read heat forcing corresponding to next twind
    if (iint == 1 .OR. mod(iint,iwater) == 0.) then
        do i=1,im
            do j=1,jm
                wssurfb(i,j)=wssurff(i,j)
            end do
        end do
        call read_water_pnetcdf((iint+iwater)/iwater,wssurff)
    end if

! linear interpolation in time
    ntime=time/twater
    fnew=time/twater-ntime
    fold=1.-fnew
    do i=1,im
        do j=1,jm
            wssurf(i,j)=fold*wssurfb(i,j)+fnew*wssurff(i,j)
        end do
    end do

    return
    end subroutine water


    subroutine fluxes1 !KTSIARAS
!!  Read and interpolate (in time) wind stress and atmospheric data.
    implicit none
    include 'pom.h'
    integer :: nfre,nstart
    parameter (nfre=3,nstart=6) !3 hour, 18UTC
    integer :: i,j,ntime,iflux,n,iold,inew
!      real uair(im,jm),vair(im,jm),tair(im,jm),rhum(im,jm),
!     & sol(im,jm),rlond(im,jm),rain(im,jm)

    real :: tflux,fold,fnew,tfluxd

    tflux=nfre ! time between wind files (hours)
    tfluxd=(nfre*3600./86400.) !(days)
    iflux=int(tflux*3600.e0/dti)

! read wind stress data
! read initial wind file
    if (iint == 1) then
        call read_atmos_pnetcdf(iint/iflux+1,uairf,vairf, &
        tairf,rhumf,solf,rlondf,rainf) !KTSIARAS
    endif

    if (iint == 1)then
        do n=1,9
            write(*,*)'TEST-TAIR',time,my_task,rhumf(29,40,n)
        enddo
    endif
! read wind file corresponding to next time
    if (mod(iint,iflux) == 0.) iold=iold+1

! linear interpolation in time
    ntime=(time-int(time0))/tfluxd
    iold=ntime-nstart+1
    inew=iold+1
    if(inew > 9)inew=9
    fnew=(time-int(time0))/tfluxd-ntime
    fold=1.-fnew

    do i=1,im
        do j=1,jm
            uair(i,j)=fold*uairf(i,j,iold)+fnew*uairf(i,j,inew)
            vair(i,j)=fold*vairf(i,j,iold)+fnew*vairf(i,j,inew)
            tair(i,j)=fold*tairf(i,j,iold)+fnew*tairf(i,j,inew)
            rhum(i,j)=fold*rhumf(i,j,iold)+fnew*rhumf(i,j,inew)
            sol(i,j)=fold*solf(i,j,iold)+fnew*solf(i,j,inew)
            rlond(i,j)=fold*rlondf(i,j,iold)+fnew*rlondf(i,j,inew)
            rain(i,j)=fold*rainf(i,j,iold)+fnew*rainf(i,j,inew)
        end do
    end do
    write(*,*)'TEST-INTERP',time,my_task,ntime,tfluxd,fnew,fold,inew, &
    iold,rhum(29,40),rhumf(29,40,1),rhumf(29,40,2)
!--calculate fluxes (wusurf,wvsurf,wtsurf,emp) from bulk-formulas
    if(iint == iend)then
        CALL PRXY('rhum',TIME,rhum(1,1),IM,1,JM,1,1.E0)
    endif

!      call bulk

    return
    end subroutine fluxes1



    subroutine runoff !KTSIARAS
    include 'pom.h'
    include 'rivers.h'
    integer :: iriv,mp12,mm12
    real :: xp12,xm12
    real :: qyr(nriver)
    data rper,rshift/365.,-60./

    CALL CLOCK12(time,2000.,mp12,mm12,xp12,xm12)

!--calculate freshwater input
!-- roff,rivqmo,ano_q_rivT specified in rivers.h
    do j=1,jm
        do i=1,im
            iriv=int(riv(i,j))
            if(iriv /= 0 .AND. iriv /= 99 .AND. iriv /= 199)then
                nr=iriv
            !        write(*,*)'TEST-RIV',my_task,i,j,nr,npoi(nr)
                qyr(nr)=roff(nr) !annual runoff
                if(rivqmo(nr,mm12) > 0.)then
                    qval(nr) = qyr(nr)*(rivqmo(nr,mm12)*xm12 + &
                    rivqmo(nr,mp12)*xp12)
                else
                    qval(nr) = qyr(nr)+.6*roff(nr)*cos(2*pi*(rshift+time)/rper)
                endif
                div=1./float(npoi(nr))
                rr=div*qval(nr)
            !        write(*,*)'TEST-Q',time,qval(nr),div,art(i,j)
                vfluxf(i,j)= vfluxf(i,j)  -rr/art(i,j)
                vfluxb(i,j)= vfluxf(i,j)
            endif
        enddo
    enddo
!      write(*,*)'OK-RIV'
    return
    end subroutine runoff

    subroutine dard(idx) !KTSIARAS
!!  Calculate water input, temperature, salinity for Dardanelles open boundary condition.
    include 'pom.h'
    include 'rivers.h'
    real :: qyr(nriver)
    data rper,rshift/365.,-60./
    integer :: iflag,idx,i,j,ii,jj
    DATA SBIAS/35.0/
    DATA TBIAS/10.0/
         
    CALL CLOCK12(time,2000.,mp12,mm12,xp12,xm12)

!--calculate freshwater input
!-- roff,rivqmo,ano_q_rivT specified in rivers.h
    iflag=0
    do j=1,jm
        do i=1,im
            iriv=int(riv(i,j))
            if(iriv == 99)then
                ii=i
                jj=j
                iflag=1
            endif
        enddo
    enddo

    if(iflag == 0)return

    i=ii
    j=jj
    QIN= 0.032*1000000.
    QOUT=0.022*1000000.
    QPLA=0.005*1000000.
! --- Modify qin and qout (seasonal cycle)
    DARDPER=365.*24*3600.
    TSHIF1=-75.*24.*3600.
    TSHIF2= -120.*24.*3600.
    TAF1=2*PI*(TSHIF1+TIME*24.*3600.)
    TAF2=2*PI*(TSHIF2+TIME*24.*3600.)
    QIN =QIN +QPLA*COS(TAF1/DARDPER)
    QOUT=QOUT+QPLA*SIN(TAF2/DARDPER)
! ------ Salinity and Temperature Seasonal Cycle -----
    SPLA=2.0
    TSHIF3=-30.*24.*3600.
    TAF3=2*PI*(TSHIF3+TIME*24.*3600.)
    sdard=28.5+SPLA*COS(TAF3/DARDPER)
    sdard=sdard-SBIAS
! old sdard=29.6-SBIAS
    saeg=38.8-SBIAS
    TPLA=6.0
    TSHIF4=-160.*24.*3600.
    TAF4=2*PI*(TSHIF4+TIME*24.*3600.)
!      tdard=16.0+TPLA*SIN(TAF4/DARDPER)
!--tdard using monthly values from Besiktepe 1993 !OLD
!--tdard using AVHRR sst
    tdard = dard_sst(mm12)*xm12 + dard_sst(mp12)*xp12
    tdard=tdard-TBIAS

    DEP1=25.
    DEP2=H(i,j)-DEP1

    EMV1=DY(i,j)*DEP1
    EMV2=DY(i,j)*DEP2

    UIN=QIN/EMV1
    UOUT=QOUT/EMV2
    UNET=(QIN-QOUT)/(EMV1+EMV2)

    GO TO (10,20,30,40,50,60) IDX

    10 CONTINUE
    ELF(I,J)=ELF(I-1,J)
    GOTO 1000

! --- EXTERNAL VELOCITY FOR RIVER
    20 CONTINUE
    UAF(I,J)=-UNET
    VAF(I,J)=0.0

    GOTO 1000

! --- VELOCITY FOR RIVER
    30 CONTINUE
    DO K=1,KBM1
        dep=-ZZ(k)*DT(i,j)
        if(dep <= dep1) THEN
            UF(I,J,K)=-UIN
        ELSE
            UF(I,J,K)=UOUT
        END IF
        VF(I,J,K)=0.0
    END DO
    GOTO 1000
! --- SALINITY & TEMPERATURE FOR RIVER
    40 CONTINUE
    DO K=1,KBM1
        DEP=-ZZ(k)*DT(I,J)
        IF(DEP <= DEP1) THEN
            UF(I,J,K)=tdard
            VF(I,J,K)=sdard
        ELSE
            VF(I,J,K)=saeg
            UF(I,J,K)=UF(I-1,J,K)
        END IF
    !         UF(I,J,K)=UF(I-1,J,K)  ! Temperature
    END DO
    if(iint == 1)write(*,*)'TEST-TDARD',tdard+10., &
    dard_sst(mm12),dard_sst(mp12),xm12,xp12,UF(I,J,1)

    GOTO 1000

! --- VERTICAL VELOCITY FOR RIVER
    50 CONTINUE
    DO K=1,KB
        W(I,J,K)=0.E0
    END DO
    GOTO 1000

! --- Q2 & L FOR RIVER
    60 CONTINUE
    DO K=1,KB
        UF(I,J,K)=1.E-10
        VF(I,J,K)=1.E-10
    END DO
    GOTO 1000


    1000 CONTINUE

    110 CONTINUE

    RETURN
    end subroutine dard


    subroutine bulk !KTSIARAS
!!  SUBROUTINE TO CALCULATE SURFACE FLUXES FROM ETA ATMOS. MODEL DATA  
!!  WUSURF, WVSURF, WTSURF (LOSS), EMP ARE CALCULATED IN BULK         
!!  SWRAD IS GIVEN DIRECTLY FROM ETA MODEL                             
!!  DOWNWARD LONGWAVE RADIATION IS GIVEN DIRECTLY FROM ETA MODEL      
!!  UPWARD LONGWAVE RADIATION IS GIVEN BY BIGNAMI 1995
!!  UNITS IN MKS                                                     
!!  FORMULAS TAKEN FROM CASTELARI ET.AL. 1992                       
    INCLUDE 'pom.h'
    real :: SREL(IM,JM),gzs(im,jm)

!--------------------------------------------------------------------
!       coefficients ( in MKS )  :
!-----------------------------------------------------------------

! --- surface air pressure, expsi, dry air gas constant

    data ps,expsi,rd / 1013., 0.622, 287./

! --- turbulent exchange coefficients ( from Rosati & Miyakoda 1988 )

    data  ce1,ch1  / 1.1e-3, 1.1e-3/

! --- turbulent exchange coefficients ( from Budyko 1963 )

    data  ce2,ch2  / 2.1e-3, 2.1e-3/

! --- air density, Stefan-Boltzmann constant , ocean emissivity

    data arho,sigma ,emiss   /1.2,   5.67e-8, .97/

! --- Solar constant , specific heat capacity

    data solar ,cp   /1350., 1005./


    data ckelv /273.16/

    const = 0.622/1013.
    tday=3600.*24.
    gz=2.0/tday

! --- Define a space varying gz for S and T --------------
    do  i=1,im
        do  j=1,jm
            gzs(i,j)=.20/tday
        enddo
    enddo

    DO 1000 i=1,im
        DO 1000 j=1,jm
            if(fsm(i,j) == 0.)cycle
            unow = uair(i,j)
            vnow = vair(i,j)
            tnow = tair(i,j)
            rhnow = rhum(i,j)
            precip= rain(i,j)
            sst =   t(i,j,1)+10.
            dlong = rlond(i,j)
            sola = sol(i,j)

            sst =   t(i,j,1)+10.
            sss =   s(i,j,1)+35.
            sss0=   sclim(i,j,1)+35.
        
        
        
        ! --- compute sp for bulk coeff.
        
            SP = sqrt(unow*unow+vnow*vnow)

        !---SP is already estimated in FORCING.ECMWF------

        
        ! --- SST data converted in Kelvin degrees
        ! --- TNOW is already in Kelvin
        
            sstk = sst + ckelv
            tnowk = tnow
        
        
        ! ---calculates the Saturation Vapor Pressure at air temp. and at sea temp.
        ! ---esat(Ta) , esat(Ts)
        
            esatair = esk(tnowk)
            esatoce = esk(sstk)
        
        ! --- calculates the saturation mixing ratios at air temp. and sea temp.
        ! --- wsat(Ta) , wsat(Ts)
        
            wsatair = (expsi/ps) * esatair
            wsatoce = (expsi/ps) * esatoce
        
        ! --- calculates the mixing ratio of the air
        ! --- w(Ta)
        
            wair = 0.01 * rhnow * wsatair
        
        ! --- calculates the density of the moist air
        
            rhom = 100.*(ps/rd) * (expsi*(1.+wair)/(tnowk*(expsi+wair)))
        
        !---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        !       QBW  (the long wave flux)
        !---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        
        
        ! --- Bigniami Original formula
        
        !            ea12 = sqrt(0.01*rhnow*esatair)
        
        !        QBW=0.98*sigma*sstk**4*(0.344-6.66e-3*ea12)*(1.-0.42*cld)


        !  ----- Bignami et al. 1995
        
            QBW=0.98*sigma*sstk**4 - dlong

        
        ! --- May formula from May (1986) :
        
        !            ea12 = sqrt(0.01*rhnow*esatair)

        !            QBW = (1.-0.75*(cld**3.4))
        !    &             * (sigma*(tnowk**4.)*(0.4 -0.05*ea12)
        !    &             + 4.*sigma*(tnowk**3.)*(sstk-tnowk))

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        ! --- calculates the term : ( Ts - Ta )
            deltemp = sstk - tnowk
        
        !---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        !       QH  (the sensible heat flux) and  QE  (the latent heat flux)
        !---- ----- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        
        ! --- calculates the QH and QE with the constant values for Ch,Ce from
        ! --- Budyko (1963)
        
        !           QH = rhom*cp*ch2*sp*deltemp
        ! --- calculates the QH and QE with variable values for Ch,Ce from
        ! --- Kondo (1975)
        
        ! --- calculate S :
        
            ss=0.
            if(sp > 0.0)ss=deltemp/(sp**2.)
        
        ! --- calculate the Stability Parameter :
        
            stp=ss*abs(ss)/(abs(ss)+0.01)
        
        ! --- for stable condition :
        
            IF(ss < 0.) THEN
                if((stp > -3.3) .AND. (stp < 0.)) then
                    fh=0.1+0.03*stp+0.9*exp(4.8*stp)
                    fe=fh
                else
                    if(stp <= -3.3) then
                        fh=0.
                        fe=fh
                    
                    endif
                endif
            
            ! --- for unstable condition :
            
            ELSE
                fh=1.0+0.63*sqrt(stp)
                fe=fh
            ENDIF
        
        ! --- calculate the coefficient CH3,CE3,CD3
        
            ch3=1.3e-03*fh
            ce3=1.5e-03*fe
        
        
            QH = rhom*ch3*cp*sp*deltemp
        
        
        ! --- calculates the term : esat(Ts)-r*esat(Ta)
        
            evap  = esatoce - rhnow*0.01*esatair
        
        ! --- calculates the term : Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/1013
        ! --- Evaporation rate [kg/(m2*sec)]
        
            EVAP = rhom*ce3*sp*evap*const
        
        ! --- calculates the LATENT HEAT FLUX  QE in MKS ( watt/m*m )
        ! --- QE = L*rho*Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/1013
        
            QE = evap*heatlat(sst)
            if(i == 340. .AND. j == 40)then
                write(*,*)'FLU',qe,evap,sst,heatlat(sst),ce3,sp,rhnow
            endif

        
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        !---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        !       WFLUX  ( the water flux ) in m/sec
        !---- -- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        
            WFLUX =  evap/1023. - precip
        
        ! ---- here relax surface salinity to seasonal climatology
        ! %%%%%%%%% RELAXATION OF SALINITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! %%%%%%%%% REMEMBER TO CHANGE THE FORMULA BEFORE PROFT %%%%%%%%%%%%
            SFLUX = (sss0-sss)*gzs(i,j)
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        ! --- Reverse Sign for model needs
        ! ifdef freshwater
        !              emp(i,j)= wflux + sflux/(sss+35.)
            vfluxf(i,j)= wflux + sflux/(sss+35.)
            vfluxb(i,j)= vfluxf(i,j)
            empbio(i,j)= -vfluxf(i,j)

            wssurf(i,j)=0.
        ! else
        !              emp(i,j)= -wflux
        !              srel(i,j)=-sflux
        ! endif
        
        !---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        !       QU ( the net upward flux )
        !---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        
        ! --- calculates : Qu = Qb + QH + QE  ( Rosati,Miyakoda 1988 ; eq. 3.2 )
        
            QU = qbw + qh + qe
        
        ! --- Divide by 1023*3991 (rho*cp) for model needs
        ! --- DO NOT Reverse Sign because it is already positive
        
            wtsurf(i,j) = QU*(1./4.082793e6)
            swrad(i,j) = -sola/4.082793e6
        !              wtsurf(i,j) =0.
        !              swrad(i,j) = 0.
                          
        !---- ----- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        !       TAUX, TAUY ( the wind stresses )
        !---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        
        ! --- calculates the Drag Coefficient Cd with variable value from
        ! --- Hellerman & Rosenstein (1983)
        
        ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! --- (5-10-92) we discover that CD of H. R. needs to be used with
        ! --- deltem = Ts - Ta
        ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            cd1 = cd(sp,deltemp)
        
        ! --- calculates the wind stresses in MKS ( newton/m*m )
        ! --- taux= rho*Cd*|V|u     tauy= rho*Cd*|V|v
        
            TAUX= rhom*cd1*sp*unow
            TAUY= rhom*cd1*sp*vnow
        ! --- Reverse Sign and devide by 1023 (rho) for model needs
            wusurf(i,j)=-taux/1023.
            wvsurf(i,j)=-tauy/1023.
        !             wusurf(i,j)=0.
        !             wvsurf(i,j)=0.

        !             if(i.eq.1.and.j.eq.54)write(*,*)'TEST-NAN',wusurf(i,j),
        !     & TAUX,rhom,cd1,sp,unow

        
        ! ------- MULTIPLY BY MASK  -------------------------------------
        
            wusurf(i,j)=wusurf(i,j)*fsm(i,j)
            wvsurf(i,j)=wvsurf(i,j)*fsm(i,j)
            wtsurf(i,j)=wtsurf(i,j)*fsm(i,j)
            vfluxf(i,j)=vfluxf(i,j)*fsm(i,j)
            vfluxb(i,j)=vfluxb(i,j)*fsm(i,j)

        !      emp(i,j)=emp(i,j)*fsm(i,j)
            srel(i,j)=srel(i,j)*fsm(i,j)


    1000 END DO

!      if(iint.lt.5)then
!      CALL PRXY('UAIR',TIME,uair(1,1),IM,1,JM,1,1.E0)
!      CALL PRXY('TAIR',TIME,tair(1,1),IM,1,JM,1,1.E0)
!      CALL PRXY('RHUM',TIME,rhum(1,1),IM,1,JM,1,1.E0)
!      CALL PRXY('SOL',TIME,sol(1,1),IM,1,JM,1,1.E0)
!      endif

!      write(*,*)'FLUX',wtsurf(232,88),swrad(232,88),wusurf(232,88)
!      write(*,*)'ATMOS',fsm(232,88),sol(232,88),tair(232,88),
!     &uair(232,88)

! --------------------------------------------------
    RETURN

    end subroutine bulk

    function HEATLAT(t) !KTSIARAS

!!  Calculates the Latent Heat of Vaporization ( J/kg ) as function of
!!  the temperature ( Celsius degrees )
!!  ( from A. Gill  pag. 607 )
! --- Constant Latent Heat of Vaporization ( Rosati,Miyakoda 1988 )
!     L = 2.501e+6  (MKS)

    heatlat = 2.5008e+6 -2.3e+3*t

    return
    end function HEATLAT

    function CD(sp,delt)

!!  Calculates the Drag Coefficient as a function of the abs. value of the
!!  wind velocity
!!  ( Hellermann and  Rosenstein )
    dimension a(6)
    data a / 0.934e-3,0.788e-4,0.868e-4 &
    ,-0.616e-6,-.120e-5,-.214e-5/

    cd = a(1) + a(2)*sp + a(3)*delt + a(4)*sp*sp &
    + a(5)*delt*delt  + a(6)*sp*delt

    return
    end function CD

!==============================================================================

    real function ESK(t) !KTSIARAS

!! Compute the saturation water vapor pressure from
!! temperature in kelvin  (Lowe,1977; JAM,16,100,1977)
    dimension a(7)
    data  a /6984.505294,-188.9039310,2.133357675,-1.288580973e-2, &
    4.393587233e-5,-8.023923082e-8,6.136820929e-11  /

    esk = a(1) +t*(a(2)+t*(a(3)+t*(a(4)+t*(a(5)+t*(a(6)+t*a(7))))))

    return
    end function ESK

! _____________________________________________________________________
    subroutine restore_interior
!!  Read, interpolate (in time) and apply restore interior data
    implicit none
    include 'pom.h'
    integer :: nz
    parameter(nz=40)
    real :: z0(nz),f0(im,jm,nz)
    integer :: i,j,k,ntime,irst
    real :: trst,fold,fnew

    trst=30 ! time between restore files (days)
    irst=int(trst*86400.e0/dti)
    ntime=time/trst

! read restore data
! read initial restore file
    if (iint == 2) then
        call read_restore_t_interior_pnetcdf(iint/irst,nz,z0,f0)
        call ztosig(z0,f0,zz,h,trstrf,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        call read_restore_s_interior_pnetcdf(iint/irst,nz,z0,f0)
        call ztosig(z0,f0,zz,h,srstrf,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
        call read_restore_tau_interior_pnetcdf(iint/irst,nz,z0,f0)
        call ztosig(z0,f0,zz,h,taurstrf,im,jm,nz,kb, &
        n_west,n_east,n_south,n_north)
    end if
! read restore file corresponding to next time
    if (iint == 2 .OR. mod(iint,irst) == 0.) then
        do k=1,kbm1
            do i=1,im
                do j=1,jm
                    trstrb(i,j,k)=trstrf(i,j,k)
                    srstrb(i,j,k)=srstrf(i,j,k)
                    taurstrb(i,j,k)=taurstrf(i,j,k)
                end do
            end do
        end do
        if (iint /= iend) then
            call read_restore_t_interior_pnetcdf((iint+irst)/irst,nz,z0, &
            f0)
            call ztosig(z0,f0,zz,h,trstrf,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            call read_restore_s_interior_pnetcdf((iint+irst)/irst,nz,z0, &
            f0)
            call ztosig(z0,f0,zz,h,srstrf,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
            call read_restore_tau_interior_pnetcdf((iint+irst)/irst,nz, &
            z0,f0)
            call ztosig(z0,f0,zz,h,taurstrf,im,jm,nz,kb, &
            n_west,n_east,n_south,n_north)
        end if
    end if

! linear interpolation in time
    fnew=time/trst-ntime
    fold=1.-fnew
    do k=1,kbm1
        do i=1,im
            do j=1,jm
                trstr(i,j,k)=fold*trstrb(i,j,k)+fnew*trstrf(i,j,k)
                srstr(i,j,k)=fold*srstrb(i,j,k)+fnew*srstrf(i,j,k)
                taurstr(i,j,k)=fold*taurstrb(i,j,k)+fnew*taurstrf(i,j,k)
            end do
        end do
    end do

! restore
    do k=1,kbm1
        do i=1,im
            do j=1,jm
                t(i,j,k)=t(i,j,k)+2.*dti/86400.*taurstr(i,j,k)* &
                (trstr(i,j,k)-t(i,j,k))
                tb(i,j,k)=tb(i,j,k)+2.*dti/86400.*taurstr(i,j,k)* &
                (trstr(i,j,k)-tb(i,j,k))
                s(i,j,k)=s(i,j,k)+2.*dti/86400.*taurstr(i,j,k)* &
                (srstr(i,j,k)-s(i,j,k))
                sb(i,j,k)=sb(i,j,k)+2.*dti/86400.*taurstr(i,j,k)* &
                (srstr(i,j,k)-sb(i,j,k))
            end do
        end do
    end do

! mask
    do k=1,kbm1
        do j=1,jm
            do i=1,im
                t(i,j,k)=t(i,j,k)*fsm(i,j)
                tb(i,j,k)=tb(i,j,k)*fsm(i,j)
                s(i,j,k)=s(i,j,k)*fsm(i,j)
                sb(i,j,k)=sb(i,j,k)*fsm(i,j)
            end do
        end do
    end do

    return
    end subroutine restore_interior
