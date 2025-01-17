! solver.f

! main subroutines for POM

! ______________________________________________________________________
    subroutine advave
!! Calculate horizontal and vertical advection and diffusion. First we calculate horizontal advective fluxes, then we add viscous fluxes.
    implicit none
    include 'pom.h'
    real :: curv2d(im,jm)
    integer :: i,j

! u-advection and diffusion

! advective fluxes
    do j=1,jm
        do i=1,im
            advua(i,j)=0.e0
        end do
    end do

    do j=2,jm
        do i=2,imm1
            fluxua(i,j)=.125e0*((d(i+1,j)+d(i,j))*ua(i+1,j) &
            +(d(i,j)+d(i-1,j))*ua(i,j)) &
            *(ua(i+1,j)+ua(i,j))
        end do
    end do

    do j=2,jm
        do i=2,im
            fluxva(i,j)=.125e0*((d(i,j)+d(i,j-1))*va(i,j) &
            +(d(i-1,j)+d(i-1,j-1))*va(i-1,j)) &
            *(ua(i,j)+ua(i,j-1))
        end do
    end do

! add viscous fluxes
    do j=2,jm
        do i=2,imm1
            fluxua(i,j)=fluxua(i,j) &
            -d(i,j)*2.e0*aam2d(i,j)*(uab(i+1,j)-uab(i,j)) &
            /dx(i,j)
        end do
    end do

    do j=2,jm
        do i=2,im
            tps(i,j)=.25e0*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1)) &
            *(aam2d(i,j)+aam2d(i,j-1) &
            +aam2d(i-1,j)+aam2d(i-1,j-1)) &
            *((uab(i,j)-uab(i,j-1)) &
            /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1)) &
            +(vab(i,j)-vab(i-1,j)) &
            /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
            fluxua(i,j)=fluxua(i,j)*dy(i,j)
            fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.25e0 &
            *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
    end do

    call exchange2d_mpi(fluxua,im,jm)
!      call exchange2d_mpi(fluxva,im,jm) !KTSIARAS

    do j=2,jmm1
        do i=2,imm1
            advua(i,j)=fluxua(i,j)-fluxua(i-1,j) &
            +fluxva(i,j+1)-fluxva(i,j)
        end do
    end do

! v-advection and diffusion
    do j=1,jm
        do i=1,im
            advva(i,j)=0.e0
        end do
    end do

! advective fluxes
    do j=2,jm
        do i=2,im
            fluxua(i,j)=.125e0*((d(i,j)+d(i-1,j))*ua(i,j) &
            +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1)) &
            *(va(i-1,j)+va(i,j))
        end do
    end do

    do j=2,jmm1
        do i=2,im
            fluxva(i,j)=.125e0*((d(i,j+1)+d(i,j))*va(i,j+1) &
            +(d(i,j)+d(i,j-1))*va(i,j)) &
            *(va(i,j+1)+va(i,j))
        end do
    end do

! add viscous fluxes
    do j=2,jmm1
        do i=2,im
            fluxva(i,j)=fluxva(i,j) &
            -d(i,j)*2.e0*aam2d(i,j)*(vab(i,j+1)-vab(i,j)) &
            /dy(i,j)
        end do
    end do

    do j=2,jm
        do i=2,im
            fluxva(i,j)=fluxva(i,j)*dx(i,j)
            fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.25e0 &
            *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
    end do

!      call exchange2d_mpi(fluxua,im,jm) !KTSIARAS
    call exchange2d_mpi(fluxva,im,jm)

    do j=2,jmm1
        do i=2,imm1
            advva(i,j)=fluxua(i+1,j)-fluxua(i,j) &
            +fluxva(i,j)-fluxva(i,j-1)
        end do
    end do

!      call exchange2d_mpi(advua,im,jm) !KTSIARAS
!      call exchange2d_mpi(advva,im,jm) !KTSIARAS

    if(mode == 2) then

        do j=2,jmm1
            do i=2,imm1
                wubot(i,j)=-0.5e0*(cbc(i,j)+cbc(i-1,j)) &
                *sqrt(uab(i,j)**2 &
                +(.25e0*(vab(i,j)+vab(i,j+1) &
                +vab(i-1,j)+vab(i-1,j+1)))**2) &
                *uab(i,j)
            end do
        end do

        do j=2,jmm1
            do i=2,imm1
                wvbot(i,j)=-0.5e0*(cbc(i,j)+cbc(i,j-1)) &
                *sqrt(vab(i,j)**2 &
                +(.25e0*(uab(i,j)+uab(i+1,j) &
                +uab(i,j-1)+uab(i+1,j-1)))**2) &
                *vab(i,j)
            end do
        end do

        do j=2,jmm1
            do i=2,imm1
                curv2d(i,j)=.25e0 &
                *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j)) &
                -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1))) &
                /(dx(i,j)*dy(i,j))
            end do
        end do
        call exchange2d_mpi(curv2d,im,jm)

        do j=2,jmm1
            if(n_west == -1) then
                do i=3,imm1
                    advua(i,j)=advua(i,j)-aru(i,j)*.25e0 &
                    *(curv2d(i,j)*d(i,j) &
                    *(va(i,j+1)+va(i,j)) &
                    +curv2d(i-1,j)*d(i-1,j) &
                    *(va(i-1,j+1)+va(i-1,j)))
                end do
            else
                do i=2,imm1
                    advua(i,j)=advua(i,j)-aru(i,j)*.25e0 &
                    *(curv2d(i,j)*d(i,j) &
                    *(va(i,j+1)+va(i,j)) &
                    +curv2d(i-1,j)*d(i-1,j) &
                    *(va(i-1,j+1)+va(i-1,j)))
                end do
            end if
        end do

        do i=2,imm1
            if(n_south == -1) then
                do j=3,jmm1
                    advva(i,j)=advva(i,j)+arv(i,j)*.25e0 &
                    *(curv2d(i,j)*d(i,j) &
                    *(ua(i+1,j)+ua(i,j)) &
                    +curv2d(i,j-1)*d(i,j-1) &
                    *(ua(i+1,j-1)+ua(i,j-1)))
                end do
            else
                do j=2,jmm1
                    advva(i,j)=advva(i,j)+arv(i,j)*.25e0 &
                    *(curv2d(i,j)*d(i,j) &
                    *(ua(i+1,j)+ua(i,j)) &
                    +curv2d(i,j-1)*d(i,j-1) &
                    *(ua(i+1,j-1)+ua(i,j-1)))
                end do
            end if
        end do

    end if

    return
    end subroutine advave

! ______________________________________________________________________
    subroutine advct
!! Calculate the horizontal portions of momentum advection well in
!! advance of their use in advu and advv so that their vertical integrals
!! (created in the main program) may be used in the external (2-D) mode
!! calculation
    implicit none
    include 'pom.h'
    real :: xflux(im,jm,kb),yflux(im,jm,kb)
    real :: curv(im,jm,kb)
    real :: dtaam
    integer :: i,j,k

    do k=1,kb
        do j=1,jm
            do i=1,im
                curv(i,j,k)=0.e0
                advx(i,j,k)=0.e0
                xflux(i,j,k)=0.e0
                yflux(i,j,k)=0.e0
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                curv(i,j,k)=.25e0*((v(i,j+1,k)+v(i,j,k)) &
                *(dy(i+1,j)-dy(i-1,j)) &
                -(u(i+1,j,k)+u(i,j,k)) &
                *(dx(i,j+1)-dx(i,j-1))) &
                /(dx(i,j)*dy(i,j))
            end do
        end do
    end do
    call exchange3d_mpi(curv(:,:,1:kbm1),im,jm,kbm1)

! calculate x-component of velocity advection

! calculate horizontal advective fluxes
    do k=1,kbm1
        do j=1,jm
            do i=2,imm1
                xflux(i,j,k)=.125e0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k) &
                +(dt(i,j)+dt(i-1,j))*u(i,j,k)) &
                *(u(i+1,j,k)+u(i,j,k))
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jm
            do i=2,im
                yflux(i,j,k)=.125e0*((dt(i,j)+dt(i,j-1))*v(i,j,k) &
                +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k)) &
                *(u(i,j,k)+u(i,j-1,k))
            end do
        end do
    end do

! add horizontal diffusive fluxes
    do k=1,kbm1
        do j=2,jm
            do i=2,imm1
                xflux(i,j,k)=xflux(i,j,k) &
                -dt(i,j)*aam(i,j,k)*2.e0 &
                *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
                dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1)) &
                *(aam(i,j,k)+aam(i-1,j,k) &
                +aam(i,j-1,k)+aam(i-1,j-1,k))
                yflux(i,j,k)=yflux(i,j,k) &
                -dtaam*((ub(i,j,k)-ub(i,j-1,k)) &
                /(dy(i,j)+dy(i-1,j) &
                +dy(i,j-1)+dy(i-1,j-1)) &
                +(vb(i,j,k)-vb(i-1,j,k)) &
                /(dx(i,j)+dx(i-1,j) &
                +dx(i,j-1)+dx(i-1,j-1)))

                xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
                yflux(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j) &
                +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
            end do
        end do
    end do

    call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)
!      call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1) !KTSIARAS

! do horizontal advection
    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k) &
                +yflux(i,j+1,k)-yflux(i,j,k)
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            if(n_west == -1) then
                do i=3,imm1
                    advx(i,j,k)=advx(i,j,k) &
                    -aru(i,j)*.25e0 &
                    *(curv(i,j,k)*dt(i,j) &
                    *(v(i,j+1,k)+v(i,j,k)) &
                    +curv(i-1,j,k)*dt(i-1,j) &
                    *(v(i-1,j+1,k)+v(i-1,j,k)))
                end do
            else
                do i=2,imm1
                    advx(i,j,k)=advx(i,j,k) &
                    -aru(i,j)*.25e0 &
                    *(curv(i,j,k)*dt(i,j) &
                    *(v(i,j+1,k)+v(i,j,k)) &
                    +curv(i-1,j,k)*dt(i-1,j) &
                    *(v(i-1,j+1,k)+v(i-1,j,k)))
                end do
            end if
        end do
    end do

! calculate y-component of velocity advection

    do k=1,kb
        do j=1,jm
            do i=1,im
                advy(i,j,k)=0.e0
                xflux(i,j,k)=0.e0
                yflux(i,j,k)=0.e0
            end do
        end do
    end do

! calculate horizontal advective fluxes
    do k=1,kbm1
        do j=2,jm
            do i=2,im
                xflux(i,j,k)=.125e0*((dt(i,j)+dt(i-1,j))*u(i,j,k) &
                +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k)) &
                *(v(i,j,k)+v(i-1,j,k))
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=1,im
                yflux(i,j,k)=.125e0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k) &
                +(dt(i,j)+dt(i,j-1))*v(i,j,k)) &
                *(v(i,j+1,k)+v(i,j,k))
            end do
        end do
    end do

! add horizontal diffusive fluxes
    do k=1,kbm1
        do j=2,jmm1
            do i=2,im
                dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1)) &
                *(aam(i,j,k)+aam(i-1,j,k) &
                +aam(i,j-1,k)+aam(i-1,j-1,k))
                xflux(i,j,k)=xflux(i,j,k) &
                -dtaam*((ub(i,j,k)-ub(i,j-1,k)) &
                /(dy(i,j)+dy(i-1,j) &
                +dy(i,j-1)+dy(i-1,j-1)) &
                +(vb(i,j,k)-vb(i-1,j,k)) &
                /(dx(i,j)+dx(i-1,j) &
                +dx(i,j-1)+dx(i-1,j-1)))
                yflux(i,j,k)=yflux(i,j,k) &
                -dt(i,j)*aam(i,j,k)*2.e0 &
                *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)

                xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j) &
                +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
                yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
            end do
        end do
    end do

    call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)
!      call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1) !KTSIARAS

! do horizontal advection
    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k) &
                +yflux(i,j,k)-yflux(i,j-1,k)
            end do
        end do
    end do

    do k=1,kbm1
        do i=2,imm1
            if(n_south == -1) then
                do j=3,jmm1
                    advy(i,j,k)=advy(i,j,k) &
                    +arv(i,j)*.25e0 &
                    *(curv(i,j,k)*dt(i,j) &
                    *(u(i+1,j,k)+u(i,j,k)) &
                    +curv(i,j-1,k)*dt(i,j-1) &
                    *(u(i+1,j-1,k)+u(i,j-1,k)))
                end do
            else
                do j=2,jmm1
                    advy(i,j,k)=advy(i,j,k) &
                    +arv(i,j)*.25e0 &
                    *(curv(i,j,k)*dt(i,j) &
                    *(u(i+1,j,k)+u(i,j,k)) &
                    +curv(i,j-1,k)*dt(i,j-1) &
                    *(u(i+1,j-1,k)+u(i,j-1,k)))
                end do
            end if
        end do
    end do

    return
    end subroutine advct

! ______________________________________________________________________
    subroutine advq(qb,q,qf)
!!  Calculates horizontal advection and diffusion, and vertical advection
!!  for turbulent quantities
    implicit none
    include 'pom.h'
    real :: qb(im_local,jm_local,kb),q(im_local,jm_local,kb)
    real :: qf(im_local,jm_local,kb)
    real :: xflux(im,jm,kb),yflux(im,jm,kb)
    integer :: i,j,k

    do k=1,kb
        do j=1,jm
            do i=1,im
                xflux(i,j,k)=0.e0
                yflux(i,j,k)=0.e0
            end do
        end do
    end do

! do horizontal advection
    do k=2,kbm1
        do j=2,jm
            do i=2,im
                xflux(i,j,k)=.125e0*(q(i,j,k)+q(i-1,j,k)) &
                *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
                yflux(i,j,k)=.125e0*(q(i,j,k)+q(i,j-1,k)) &
                *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
            end do
        end do
    end do

! do horizontal diffusion
    do k=2,kbm1
        do j=2,jm
            do i=2,im
                xflux(i,j,k)=xflux(i,j,k) &
                -.25e0*(aam(i,j,k)+aam(i-1,j,k) &
                +aam(i,j,k-1)+aam(i-1,j,k-1)) &
                *(h(i,j)+h(i-1,j)) &
                *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j) &
                /(dx(i,j)+dx(i-1,j))
                yflux(i,j,k)=yflux(i,j,k) &
                -.25e0*(aam(i,j,k)+aam(i,j-1,k) &
                +aam(i,j,k-1)+aam(i,j-1,k-1)) &
                *(h(i,j)+h(i,j-1)) &
                *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j) &
                /(dy(i,j)+dy(i,j-1))
                xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
                yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
            end do
        end do
    end do

    call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)
    call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

! do vertical advection, add flux terms, then step forward in time
    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1)) &
                *art(i,j)/(dz(k)+dz(k-1)) &
                +xflux(i+1,j,k)-xflux(i,j,k) &
                +yflux(i,j+1,k)-yflux(i,j,k)
                qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j) &
                *qb(i,j,k)-dti2*qf(i,j,k)) &
                /((h(i,j)+etf(i,j))*art(i,j))
            end do
        end do
    end do

    return
    end subroutine advq

! ______________________________________________________________________
    subroutine advt1(fb,f,fclim,ff)
!!  Integrate conservative scalar equations.
!!  This is using a centred scheme, as originally provided in POM (previously
!!  called advt)
    implicit none
    include 'pom.h'
    real :: fb(im_local,jm_local,kb),f(im_local,jm_local,kb)
    real :: fclim(im_local,jm_local,kb),ff(im_local,jm_local,kb)
    real :: xflux(im,jm,kb),yflux(im,jm,kb)
    integer :: i,j,k

    do k=1,kb
        do j=1,jm
            do i=1,im
                xflux(i,j,k)=0.e0
                yflux(i,j,k)=0.e0
            end do
        end do
    end do

    do j=1,jm
        do i=1,im
            f(i,j,kb)=f(i,j,kbm1)
            fb(i,j,kb)=fb(i,j,kbm1)
        end do
    end do

!      call exchange3d_mpi(f(:,:,1:kbm1),im,jm,kbm1) !KTSIARAS

! do advective fluxes
    do k=1,kbm1
        do j=2,jm
            do i=2,im
                xflux(i,j,k)=.25e0*((dt(i,j)+dt(i-1,j)) &
                *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
                yflux(i,j,k)=.25e0*((dt(i,j)+dt(i,j-1)) &
                *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
            end do
        end do
    end do

! add diffusive fluxes
    do k=1,kb
        do j=1,jm
            do i=1,im
                fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
            end do
        end do
    end do

!      call exchange3d_mpi(fb(:,:,1:kbm1),im,jm,kbm1) !KTSIARAS

    do k=1,kbm1
        do j=2,jm
            do i=2,im
                xflux(i,j,k)=xflux(i,j,k) &
                -.5e0*(aam(i,j,k)+aam(i-1,j,k)) &
                *(h(i,j)+h(i-1,j))*tprni &
                *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j) &
                /(dx(i,j)+dx(i-1,j))
                yflux(i,j,k)=yflux(i,j,k) &
                -.5e0*(aam(i,j,k)+aam(i,j-1,k)) &
                *(h(i,j)+h(i,j-1))*tprni &
                *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j) &
                /(dy(i,j)+dy(i,j-1))
                xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
                yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
            end do
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
            end do
        end do
    end do

!      call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1) !KTSIARAS
!      call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1) !KTSIARAS
     

! do vertical advection
    do j=2,jmm1
        do i=2,imm1
            zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
            zflux(i,j,kb)=0.e0
        end do
    end do

    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                zflux(i,j,k)=.5e0*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
            end do
        end do
    end do

! add net horizontal fluxes and then step forward in time
    do j=2,jmm1
        do i=2,imm1
            do k=1,kbm1
                ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k) &
                +yflux(i,j+1,k)-yflux(i,j,k) &
                +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
                ff(i,j,k)=(fb(i,j,k)*dble((h(i,j)+etb(i,j))*art(i,j)) &
                -dti2*ff(i,j,k))/dble((h(i,j)+etf(i,j))*art(i,j))
                if(fsm(i,j) == 1. .AND. isnan(ff(i,j,k)))then
                    write(*,*)'TESTNAN',my_task,i,j,k,u(i,j,k),v(i,j,k), &
                    elb(i,j)
                    stop
                endif
            end do
        end do
    end do

    return
    end subroutine advt1

! ______________________________________________________________________
    subroutine advt2(fb,f,fclim,ff)
!!  Integrate conservative scalar equations.
!!  This is a first-order upstream scheme, which reduces implicit
!!  diffusion using the Smolarkiewicz iterative upstream scheme with an
!!  antidiffusive velocity.
!!  It is based on the subroutines of Gianmaria Sannino (Inter-university
!!  Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
!!  National Agency for New Technology and Environment, Rome, Italy)
    implicit none
    include 'pom.h'
    real :: fb(im_local,jm_local,kb),f(im_local,jm_local,kb)
    real :: fclim(im_local,jm_local,kb),ff(im_local,jm_local,kb)
    real :: xflux(im,jm,kb),yflux(im,jm,kb)
    real :: fbmem(im,jm,kb),eta(im,jm)
    real :: xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
    integer :: i,j,k,itera

! calculate horizontal mass fluxes
    do k=1,kb
        do j=1,jm
            do i=1,im
                xflux(i,j,k)=0.e0
                yflux(i,j,k)=0.e0
                xmassflux(i,j,k)=0.e0
                ymassflux(i,j,k)=0.e0
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=2,im
                xmassflux(i,j,k)=0.25e0*(dy(i-1,j)+dy(i,j)) &
                *(dt(i-1,j)+dt(i,j))*u(i,j,k)
            end do
        end do

        do j=2,jm
            do i=2,imm1
                ymassflux(i,j,k)=0.25e0*(dx(i,j-1)+dx(i,j)) &
                *(dt(i,j-1)+dt(i,j))*v(i,j,k)
            end do
        end do
    end do

    do j=1,jm
        do i=1,im
            fb(i,j,kb)=fb(i,j,kbm1)
            eta(i,j)=etb(i,j)
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                zwflux(i,j,k)=w(i,j,k)
                fbmem(i,j,k)=fb(i,j,k)
            end do
        end do
    end do

! start Smolarkiewicz scheme
    do itera=1,nitera

    ! upwind advection scheme
        do k=1,kbm1
            do j=2,jm
                do i=2,im
                    xflux(i,j,k)=0.5e0 &
                    *((xmassflux(i,j,k)+abs(xmassflux(i,j,k))) &
                    *fbmem(i-1,j,k)+ &
                    (xmassflux(i,j,k)-abs(xmassflux(i,j,k))) &
                    *fbmem(i,j,k))

                    yflux(i,j,k)=0.5e0 &
                    *((ymassflux(i,j,k)+abs(ymassflux(i,j,k))) &
                    *fbmem(i,j-1,k)+ &
                    (ymassflux(i,j,k)-abs(ymassflux(i,j,k))) &
                    *fbmem(i,j,k))
                end do
            end do
        end do

        do j=2,jmm1
            do i=2,imm1
                zflux(i,j,1)=0.e0
                if(itera == 1) zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j)
                zflux(i,j,kb)=0.e0
            end do
        end do

        do k=2,kbm1
            do j=2,jmm1
                do i=2,imm1
                    zflux(i,j,k)=0.5e0 &
                    *((zwflux(i,j,k)+abs(zwflux(i,j,k))) &
                    *fbmem(i,j,k)+ &
                    (zwflux(i,j,k)-abs(zwflux(i,j,k))) &
                    *fbmem(i,j,k-1))
                    zflux(i,j,k)=zflux(i,j,k)*art(i,j)
                end do
            end do
        end do

    ! add net advective fluxes and step forward in time
        do j=2,jmm1
            do i=2,imm1
                do k=1,kbm1
                    ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k) &
                    +yflux(i,j+1,k)-yflux(i,j,k) &
                    +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
                    ff(i,j,k)=(fbmem(i,j,k)*dble((h(i,j)+eta(i,j))*art(i,j)) &
                    -dti2*ff(i,j,k))/dble((h(i,j)+etf(i,j))*art(i,j))
                end do
            end do
        end do
    ! next line added on 22-Jul-2009 by Raffaele Bernardello
        call exchange3d_mpi(ff(:,:,1:kbm1),im,jm,kbm1)

    ! calculate antidiffusion velocity
        call smol_adif(xmassflux,ymassflux,zwflux,ff)

        do j=1,jm
            do i=1,im
                eta(i,j)=etf(i,j)
                do k=1,kb
                    fbmem(i,j,k)=ff(i,j,k)
                end do
            end do
        end do

    ! end of Smolarkiewicz scheme
    end do

! add horizontal diffusive fluxes
    do k=1,kb
        do j=1,jm
            do i=1,im
                fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jm
            do i=2,im
                xmassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i-1,j,k))
                ymassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i,j-1,k))
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jm
            do i=2,im
                xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni &
                *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j) &
                *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j))
                yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni &
                *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j) &
                *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1))
            end do
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
            end do
        end do
    end do

! add net horizontal fluxes and step forward in time
    do j=2,jmm1
        do i=2,imm1
            do k=1,kbm1
                ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k) &
                +yflux(i,j+1,k)-yflux(i,j,k)) &
                /dble((h(i,j)+etf(i,j))*art(i,j))
            end do
        end do
    end do

    return
    end subroutine advt2

! ______________________________________________________________________
    subroutine advu
!!  Do horizontal and vertical advection of u-momentum, and include
!!  coriolis, surface slope and baroclinic terms
    implicit none
    include 'pom.h'
    integer :: i,j,k

! do vertical advection
    do k=1,kb
        do j=1,jm
            do i=1,im
                uf(i,j,k)=0.e0
            end do
        end do
    end do

    do k=2,kbm1
        do j=1,jm
            do i=2,im
                uf(i,j,k)=.25e0*(w(i,j,k)+w(i-1,j,k)) &
                *(u(i,j,k)+u(i,j,k-1))
            end do
        end do
    end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                uf(i,j,k)=advx(i,j,k) &
                +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k) &
                -aru(i,j)*.25e0 &
                *(cor(i,j)*dt(i,j) &
                *(v(i,j+1,k)+v(i,j,k)) &
                +cor(i-1,j)*dt(i-1,j) &
                *(v(i-1,j+1,k)+v(i-1,j,k))) &
                +grav*.125e0*(dt(i,j)+dt(i-1,j)) &
                *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j) &
                +(e_atmos(i,j)-e_atmos(i-1,j))*2.e0) &
                *(dy(i,j)+dy(i-1,j)) &
                +drhox(i,j,k)
            end do
        end do
    end do

!  step forward in time
    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j)) &
                *aru(i,j)*ub(i,j,k) &
                -2.e0*dti2*uf(i,j,k)) &
                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j)) &
                *aru(i,j))
            end do
        end do
    end do

    return
    end subroutine advu

! ______________________________________________________________________
    subroutine advv
!!  Do horizontal and vertical advection of v-momentum, and include
!!  coriolis, surface slope and baroclinic terms
    implicit none
    include 'pom.h'
    integer :: i,j,k

! do vertical advection
    do k=1,kb
        do j=1,jm
            do i=1,im
                vf(i,j,k)=0.e0
            end do
        end do
    end do

    do k=2,kbm1
        do j=2,jm
            do i=1,im
                vf(i,j,k)=.25e0*(w(i,j,k)+w(i,j-1,k)) &
                *(v(i,j,k)+v(i,j,k-1))
            end do
        end do
    end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                vf(i,j,k)=advy(i,j,k) &
                +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k) &
                +arv(i,j)*.25e0 &
                *(cor(i,j)*dt(i,j) &
                *(u(i+1,j,k)+u(i,j,k)) &
                +cor(i,j-1)*dt(i,j-1) &
                *(u(i+1,j-1,k)+u(i,j-1,k))) &
                +grav*.125e0*(dt(i,j)+dt(i,j-1)) &
                *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1) &
                +(e_atmos(i,j)-e_atmos(i,j-1))*2.e0) &
                *(dx(i,j)+dx(i,j-1)) &
                +drhoy(i,j,k)
            end do
        end do
    end do

! step forward in time
    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1)) &
                *arv(i,j)*vb(i,j,k) &
                -2.e0*dti2*vf(i,j,k)) &
                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1)) &
                *arv(i,j))
            end do
        end do
    end do

    return
    end subroutine advv

! ______________________________________________________________________
    subroutine baropg
!!  Calculate  baroclinic pressure gradient
    implicit none
    include 'pom.h'
    integer :: i,j,k

    do k=1,kb
        do j=1,jm
            do i=1,im
                rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
            end do
        end do
    end do

! calculate x-component of baroclinic pressure gradient
    do j=2,jmm1
        do i=2,imm1
            drhox(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i-1,j)) &
            *(rho(i,j,1)-rho(i-1,j,1))
        end do
    end do

    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhox(i,j,k)=drhox(i,j,k-1) &
                +grav*.25e0*(zz(k-1)-zz(k)) &
                *(dt(i,j)+dt(i-1,j)) &
                *(rho(i,j,k)-rho(i-1,j,k) &
                +rho(i,j,k-1)-rho(i-1,j,k-1)) &
                +grav*.25e0*(zz(k-1)+zz(k)) &
                *(dt(i,j)-dt(i-1,j)) &
                *(rho(i,j,k)+rho(i-1,j,k) &
                -rho(i,j,k-1)-rho(i-1,j,k-1))
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j)) &
                *drhox(i,j,k)*dum(i,j) &
                *(dy(i,j)+dy(i-1,j))
            end do
        end do
    end do

! calculate y-component of baroclinic pressure gradient
    do j=2,jmm1
        do i=2,imm1
            drhoy(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i,j-1)) &
            *(rho(i,j,1)-rho(i,j-1,1))
        end do
    end do

    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhoy(i,j,k)=drhoy(i,j,k-1) &
                +grav*.25e0*(zz(k-1)-zz(k)) &
                *(dt(i,j)+dt(i,j-1)) &
                *(rho(i,j,k)-rho(i,j-1,k) &
                +rho(i,j,k-1)-rho(i,j-1,k-1)) &
                +grav*.25e0*(zz(k-1)+zz(k)) &
                *(dt(i,j)-dt(i,j-1)) &
                *(rho(i,j,k)+rho(i,j-1,k) &
                -rho(i,j,k-1)-rho(i,j-1,k-1))
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1)) &
                *drhoy(i,j,k)*dvm(i,j) &
                *(dx(i,j)+dx(i,j-1))
            end do
        end do
    end do

    do k=1,kb
        do j=2,jmm1
            do i=2,imm1
                drhox(i,j,k)=ramp*drhox(i,j,k)
                drhoy(i,j,k)=ramp*drhoy(i,j,k)
            end do
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
            end do
        end do
    end do

    return
    end subroutine baropg

! ______________________________________________________________________
    subroutine baropg_mcc
!!  Calculate  baroclinic pressure gradient
!!  4th order correction terms, following McCalpin
    implicit none
    include 'pom.h'
    integer :: i,j,k
    real :: d4(im,jm),ddx(im,jm),drho(im,jm,kb),rhou(im,jm,kb)
    real :: rho4th(0:im,0:jm,kb),d4th(0:im,0:jm)

    do k=1,kb
        do j=1,jm
            do i=1,im
                rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
            end do
        end do
    end do

! convert a 2nd order matrices to special 4th order
! special 4th order case
    call order2d_mpi(d,d4th,im,jm)
    call order3d_mpi(rho,rho4th,im,jm,kb)

! compute terms correct to 4th order
    do i=1,im
        do j=1,jm
            ddx(i,j)=0.
            d4(i,j)=0.
        end do
    end do
    do k=1,kb
        do j=1,jm
            do i=1,im
                rhou(i,j,k)=0.
                drho(i,j,k)=0.
            end do
        end do
    end do

! compute DRHO, RHOU, DDX and D4
    do j=1,jm
        do i=2,im
            do k=1,kbm1
                drho(i,j,k)=(rho(i,j,k)-rho(i-1,j,k))*dum(i,j)
                rhou(i,j,k)=0.5*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j)
            end do
            ddx(i,j)=(d(i,j)-d(i-1,j))*dum(i,j)
            d4(i,j)=.5*(d(i,j)+d(i-1,j))*dum(i,j)
        end do
    end do

    if(n_west == -1) then
        do j=1,jm
            do i=3,imm1
                do k=1,kbm1
                    drho(i,j,k)=drho(i,j,k) - (1./24.)* &
                    (dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))- &
                    2*(rho(i,j,k)-rho(i-1,j,k))+ &
                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
                    rhou(i,j,k)=rhou(i,j,k) + (1./16.)* &
                    (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+ &
                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
                end do
                ddx(i,j)=ddx(i,j)-(1./24.)* &
                (dum(i+1,j)*(d(i+1,j)-d(i,j))- &
                2*(d(i,j)-d(i-1,j))+ &
                dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
                d4(i,j)=d4(i,j)+(1./16.)* &
                (dum(i+1,j)*(d(i,j)-d(i+1,j))+ &
                dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
            end do
        end do
    else
        do j=1,jm
            do i=2,imm1
                do k=1,kbm1
                    drho(i,j,k)=drho(i,j,k) - (1./24.)* &
                    (dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))- &
                    2*(rho(i,j,k)-rho(i-1,j,k))+ &
                    dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
                    rhou(i,j,k)=rhou(i,j,k) + (1./16.)* &
                    (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+ &
                    dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
                end do
                ddx(i,j)=ddx(i,j)-(1./24.)* &
                (dum(i+1,j)*(d(i+1,j)-d(i,j))- &
                2*(d(i,j)-d(i-1,j))+ &
                dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
                d4(i,j)=d4(i,j)+(1./16.)* &
                (dum(i+1,j)*(d(i,j)-d(i+1,j))+ &
                dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
            end do
        end do
    end if

! calculate x-component of baroclinic pressure gradient
    do j=2,jmm1
        do i=2,imm1
            drhox(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
    end do

    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhox(i,j,k)=drhox(i,j,k-1) &
                +grav*0.5e0*dzz(k-1)*d4(i,j) &
                *(drho(i,j,k-1)+drho(i,j,k)) &
                +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j) &
                *(rhou(i,j,k)-rhou(i,j,k-1))
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j)) &
                *drhox(i,j,k)*dum(i,j) &
                *(dy(i,j)+dy(i-1,j))
            end do
        end do
    end do

! compute terms correct to 4th order
    do i=1,im
        do j=1,jm
            ddx(i,j)=0.
            d4(i,j)=0.
        end do
    end do
    do k=1,kb
        do j=1,jm
            do i=1,im
                rhou(i,j,k)=0.
                drho(i,j,k)=0.
            end do
        end do
    end do

! compute DRHO, RHOU, DDX and D4
    do j=2,jm
        do i=1,im
            do k=1,kbm1
                drho(i,j,k)=(rho(i,j,k)-rho(i,j-1,k))*dvm(i,j)
                rhou(i,j,k)=.5*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j)
            end do
            ddx(i,j)=(d(i,j)-d(i,j-1))*dvm(i,j)
            d4(i,j)=.5*(d(i,j)+d(i,j-1))*dvm(i,j)
        end do
    end do

    if(n_south == -1) then
        do j=3,jmm1
            do i=1,im
                do k=1,kbm1
                    drho(i,j,k)=drho(i,j,k)-(1./24.)* &
                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))- &
                    2*(rho(i,j,k)-rho(i,j-1,k))+ &
                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
                    rhou(i,j,k)=rhou(i,j,k)+(1./16.)* &
                    (dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+ &
                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
                end do
                ddx(i,j)=ddx(i,j)-(1./24)* &
                (dvm(i,j+1)*(d(i,j+1)-d(i,j))- &
                2*(d(i,j)-d(i,j-1))+ &
                dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
                d4(i,j)=d4(i,j)+(1./16.)* &
                (dvm(i,j+1)*(d(i,j)-d(i,j+1))+ &
                dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
            end do
        end do
    else
        do j=2,jmm1
            do i=1,im
                do k=1,kbm1
                    drho(i,j,k)=drho(i,j,k)-(1./24.)* &
                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))- &
                    2*(rho(i,j,k)-rho(i,j-1,k))+ &
                    dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
                    rhou(i,j,k)=rhou(i,j,k)+(1./16.)* &
                    (dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+ &
                    dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
                end do
                ddx(i,j)=ddx(i,j)-(1./24)* &
                (dvm(i,j+1)*(d(i,j+1)-d(i,j))- &
                2*(d(i,j)-d(i,j-1))+ &
                dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
                d4(i,j)=d4(i,j)+(1./16.)* &
                (dvm(i,j+1)*(d(i,j)-d(i,j+1))+ &
                dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
            end do
        end do
    end if

! calculate y-component of baroclinic pressure gradient
    do j=2,jmm1
        do i=2,imm1
            drhoy(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
    end do

    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhoy(i,j,k)=drhoy(i,j,k-1) &
                +grav*0.5e0*dzz(k-1)*d4(i,j) &
                *(drho(i,j,k-1)+drho(i,j,k)) &
                +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j) &
                *(rhou(i,j,k)-rhou(i,j,k-1))
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1)) &
                *drhoy(i,j,k)*dvm(i,j) &
                *(dx(i,j)+dx(i,j-1))
            end do
        end do
    end do

    do k=1,kb
        do j=2,jmm1
            do i=2,imm1
                drhox(i,j,k)=ramp*drhox(i,j,k)
                drhoy(i,j,k)=ramp*drhoy(i,j,k)
            end do
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
            end do
        end do
    end do

    return
    end subroutine baropg_mcc

! ______________________________________________________________________
    subroutine density(si,ti,rhoo)
!!  Calculate (density-1000.)/rhoref.
!!  see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech., 609-611
!!  note: if pressure is not used in dens, buoyancy term (boygr) in
!!  subroutine profq must be changed (see note in subroutine profq)
    implicit none
    include 'pom.h'
    real :: si(im_local,jm_local,kb),ti(im_local,jm_local,kb)
    real :: rhoo(im_local,jm_local,kb)
    integer :: i,j,k
    real :: cr,p,rhor,sr,tr,tr2,tr3,tr4

    do k=1,kbm1
        do j=1,jm
            do i=1,im

                tr=ti(i,j,k)+tbias
                sr=si(i,j,k)+sbias
                tr2=tr*tr
                tr3=tr2*tr
                tr4=tr3*tr

            ! approximate pressure in units of bars
                p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5

                rhor=-0.157406e0+6.793952e-2*tr &
                -9.095290e-3*tr2+1.001685e-4*tr3 &
                -1.120083e-6*tr4+6.536332e-9*tr4*tr

                rhor=rhor+(0.824493e0-4.0899e-3*tr &
                +7.6438e-5*tr2-8.2467e-7*tr3 &
                +5.3875e-9*tr4)*sr &
                +(-5.72466e-3+1.0227e-4*tr &
                -1.6546e-6*tr2)*abs(sr)**1.5 &
                +4.8314e-4*sr*sr

                cr=1449.1e0+.0821e0*p+4.55e0*tr-.045e0*tr2 &
                +1.34e0*(sr-35.e0)
                rhor=rhor+1.e5*p/(cr*cr)*(1.e0-2.e0*p/(cr*cr))

                rhoo(i,j,k)=rhor/rhoref*fsm(i,j)

            end do
        end do
    end do

    return
    end subroutine density



SUBROUTINE PROFQold
!! ??
    INCLUDE 'pom.h'
    REAL :: KN
    DIMENSION GH(IM,JM,KB),SM(IM,JM,KB),SH(IM,JM,KB)
    DIMENSION PROD(IM,JM,KB),KN(IM,JM,KB),BOYGR(IM,JM,KB)
    DIMENSION DH(IM,JM),CC(IM,JM,KB),STF(IM,JM,KB)
    real :: a(im,jm,kb),c(im,jm,kb)
    real :: ee(im,jm,kb),gg(im,jm,kb)

!      EQUIVALENCE (A,SM),(C,SH),(PROD,KN),(TPS,DH)
!      EQUIVALENCE (DTEF,CC)
    DATA A1,B1,A2,B2,C1/0.92,16.6,0.74,10.1,0.08/
    DATA E1/1.8/,E2/1.33/,E3/1.0/
    DATA SQ/0.20/,SEF/1./

     
    DT2=DTI2


    TBIAS=10.
    SBIAS=35.

    DO 50 J=1,JM
        DO 50 I=1,IM
            DH(I,J)=H(I,J)+ETF(I,J)
    50 END DO
    DO 100 K=2,KBM1
        DO 100 J=1,JM
            DO 100 I=1,IM
                A(I,J,K)=-DT2*(KQ(I,J,K+1)+KQ(I,J,K)+2.*UMOL)*.5 &
                /(DZZ(K-1)*DZ(K)*DH(I,J)*DH(I,J))
                C(I,J,K)=-DT2*(KQ(I,J,K-1)+KQ(I,J,K)+2.*UMOL)*.5 &
                /(DZZ(K-1)*DZ(K-1)*DH(I,J)*DH(I,J))
    100 END DO
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!        DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B                    *
!                                                                      *
!***********************************************************************
!------  SURFACE AND BOTTOM B.C.S ------------
    CONST1=16.6**.6666667*SEF
    DO 90 J=1,JMM1
        DO 90 I=1,IMM1
            EE(I,J,1)=0.
            GG(I,J,1)=SQRT( (.5*(WUSURF(I,J)+WUSURF(I+1,J)))**2 &
            +(.5*(WVSURF(I,J)+WVSURF(I,J+1)))**2 )*CONST1
            UF(I,J,KB)=SQRT( (.5*(WUBOT(I,J)+WUBOT(I+1,J)))**2 &
            +(.5*(WVBOT(I,J)+WVBOT(I,J+1)))**2 )*CONST1
    90 END DO
!----- Calculate speed of sound squared ----------------------------
    DO 101 K=1,KBM1
        DO 101 J=1,JM
            DO 101 I=1,IM
                TP=T(I,J,K)+TBIAS
                SP=S(I,J,K)+SBIAS
            !      Calculate pressure in units of decibars
                P=-GRAV*1.025*ZZ(K)*DH(I,J)*.1
                CC(I,J,K)=1449.1+.00821*P+4.55*TP -.045*TP**2 &
                +1.34*(SP- 35.0)
                CC(I,J,K)=CC(I,J,K)/SQRT((1.-.01642*P/CC(I,J,K)) &
                *(1.-0.40*P/CC(I,J,K)**2))
    101 END DO
!----- Calculate buoyancy gradient ---------------------------------
    DO 102 K=2,KBM1
        DO 102 J=1,JM
            DO 102 I=1,IM
                Q2B(I,J,K)=ABS(Q2B(I,J,K))
                Q2LB(I,J,K)=ABS(Q2LB(I,J,K))
                BOYGR(I,J,K)=GRAV*((RHO(I,J,K-1)-RHO(I,J,K))/(DZZ(K-1)*DH(I,J))) &
                +GRAV**2*2.*1.025/(CC(I,J,K-1)**2+CC(I,J,K)**2)
    102 END DO
    DO 220 K=2,KBM1
        DO 220 J=1,JM
            DO 220 I=1,IM
                L(I,J,K)=Q2LB(I,J,K)/(Q2B(I,J,K)+SMALL)
                GH(I,J,K)=L(I,J,K)**2/(Q2B(I,J,K)+SMALL)*BOYGR(I,J,K)
                GH(I,J,K)=MIN(GH(I,J,K),.028)
    220 END DO
    DO 225 J=1,JM
        DO 225 I=1,IM
            L(I,J,1)=0.
            L(I,J,KB)=0.
            GH(I,J,1)=0.
            GH(I,J,KB)=0.
    225 END DO
!------ CALC. T.K.E. PRODUCTION -----------------------------------
    DO 120 K=2,KBM1
        DO 120 J=2,JMM1
            DO 120 I=1,IMM1
                PROD(I,J,K)=KM(I,J,K)*.25*SEF &
                *( (U(I,J,K)-U(I,J,K-1)+U(I+1,J,K)-U(I+1,J,K-1))**2 &
                +(V(I,J,K)-V(I,J,K-1)+V(I,J+1,K)-V(I,J+1,K-1))**2 ) &
                /(DZZ(K-1)*DH(I,J))**2
                PROD(I,J,K)=PROD(I,J,K)+KH(I,J,K)*BOYGR(I,J,K)
    120 END DO
    GHC=-6.0
    DO 110 K=1,KB
        DO 110 J=1,JM
            DO 110 I=1,IM
                STF(I,J,K)=1.
                IF(GH(I,J,K) < 0.) STF(I,J,K)=1.0-0.9*(GH(I,J,K)/GHC)**1.5
                IF(GH(I,J,K) < GHC) STF(I,J,K)=0.1
                DTEF(I,J,K)=Q2B(I,J,K)*SQRT(Q2B(I,J,K))/(B1*Q2LB(I,J,K)+SMALL) &
                *STF(I,J,K)
    110 END DO
    DO 140 K=2,KBM1
        DO 140 J=1,JM
            DO 140 I=1,IM
                GG(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-EE(I,J,K-1)) &
                -(2.*DT2*DTEF(I,J,K)+1.) )
                EE(I,J,K)=A(I,J,K)*GG(I,J,K)
                GG(I,J,K)=(-2.*DT2*PROD(I,J,K) &
                +C(I,J,K)*GG(I,J,K-1)-UF(I,J,K))*GG(I,J,K)
    140 END DO
    DO 150 K=1,KBM1
        KI=KB-K
        DO 150 J=1,JM
            DO 150 I=1,IM
                UF(I,J,KI)=EE(I,J,KI)*UF(I,J,KI+1)+GG(I,J,KI)
    150 END DO
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!        DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB                     *
!                                                                      *
!***********************************************************************
    DO 155 J=1,JM
        DO 155 I=1,IM
            IF(UF(I,J,2) < SMALL) UF(I,J,2)=SMALL
            EE(I,J,2)=0.
            GG(I,J,2)=-KAPPA*Z(2)*DH(I,J)*UF(I,J,2)
            VF(I,J,KB)=0.
    155 END DO
    DO 160 K=3,KBM1
        DO 160 J=1,JM
            DO 160 I=1,IM
                DTEF(I,J,K) =DTEF(I,J,K)*(1.+E2*((1./ABS(Z(K)-Z(1))+ &
                1./ABS(Z(K)-Z(KB))) *L(I,J,K)/(DH(I,J)*KAPPA))**2)
                GG(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-EE(I,J,K-1)) &
                -(DT2*DTEF(I,J,K)+1.))
                EE(I,J,K)=A(I,J,K)*GG(I,J,K)
                GG(I,J,K)=(DT2*(-PROD(I,J,K) &
                *L(I,J,K)*E1)+C(I,J,K)*GG(I,J,K-1)-VF(I,J,K))*GG(I,J,K)
    160 END DO
    DO 170 K=1,KB-2
        KI=KB-K
        DO 170 J=1,JM
            DO 170 I=1,IM
                VF(I,J,KI)=EE(I,J,KI)*VF(I,J,KI+1)+GG(I,J,KI)
    170 END DO
    DO 180 K=2,KBM1
        DO 180 J=1,JM
            DO 180 I=1,IM
                IF(UF(I,J,K) > SMALL .AND. VF(I,J,K) > SMALL) GOTO 180
                UF(I,J,K)=SMALL
                VF(I,J,K)=SMALL
    180 END DO
!***********************************************************************
!                                                                      *
!               THE FOLLOWING SECTION SOLVES FOR KM AND KH             *
!                                                                      *
!***********************************************************************
    COEF4=18.*A1*A1+9.*A1*A2
    COEF5=9.*A1*A2
! NOTE THAT SM,SH LIMIT TO INFINITY WHEN GH APPROACHES 0.0288
    DO 230 K=1,KB
        DO 230 J=1,JM
            DO 230 I=1,IM
                COEF1=A2*(1.-6.*A1/B1*STF(I,J,K))
                COEF2=3.*A2*B2/STF(I,J,K)+18.*A1*A2
                COEF3=A1*(1.-3.*C1-6.*A1/B1*STF(I,J,K))
                SH(I,J,K)=COEF1/(1.-COEF2*GH(I,J,K))
                SM(I,J,K)=COEF3+SH(I,J,K)*COEF4*GH(I,J,K)
                SM(I,J,K)=SM(I,J,K)/(1.-COEF5*GH(I,J,K))
    230 END DO

    DO 280 K=1,KB
        DO 280 J=1,JM
            DO 280 I=1,IM
                KN(I,J,K)=L(I,J,K)*SQRT(ABS(Q2(I,J,K)))
                KQ(I,J,K)=(KN(I,J,K)*.41*SH(I,J,K)+KQ(I,J,K))*.5
                KM(I,J,K)=(KN(I,J,K)*SM(I,J,K)+KM(I,J,K))*.5
                KH(I,J,K)=(KN(I,J,K)*SH(I,J,K)+KH(I,J,K))*.5
    280 END DO

    RETURN
    end SUBROUTINE PROFQold


! ______________________________________________________________________
    subroutine profq
!!  Solve for q2 (twice the turbulent kinetic energy), q2l (q2 x turbulent
!!  length scale), km (vertical kinematic viscosity) and kh (vertical
!!  kinematic diffusivity), using a simplified version of the level 2 1/2
!!  model of Mellor and Yamada (1982)
!!  In this version, the Craig-Banner sub-model whereby breaking wave tke
!!  is injected into the surface is included. However, we use an
!!  analytical solution to the near surface tke equation to solve for q2
!!  at the surface giving the same result as C-B diffusion. The new scheme
!!  is simpler and more robust than the latter scheme
    implicit none
    include 'pom.h'
    real :: a(im,jm,kb),c(im,jm,kb)
    real :: ee(im,jm,kb),gg(im,jm,kb)
    real :: sm(im,jm,kb),sh(im,jm,kb)
    real :: cc(im,jm,kb)
    real :: gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
    real :: prod(im,jm,kb)
    real :: a1,a2,b1,b2,c1
    real :: coef1,coef2,coef3,coef4,coef5
    real :: const1,e1,e2,ghc
    real :: p,sef,sp,tp
    real :: l0(im,jm)
    real :: cbcnst,surfl,shiw
    real :: utau2(im,jm)
    real :: df0,df1,df2
    integer :: i,j,k,ki

    data a1,b1,a2,b2,c1/0.92e0,16.6e0,0.74e0,10.1e0,0.08e0/
    data e1/1.8e0/,e2/1.33e0/
    data sef/1.e0/
    data cbcnst/100./surfl/2.e5/shiw/0.0/

    do j=1,jm
        do i=1,im
            dh(i,j)=h(i,j)+etf(i,j)
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                a(i,j,k)=0.e0
                c(i,j,k)=0.e0
                ee(i,j,k)=0.e0
                gg(i,j,k)=0.e0
            end do
        end do
    end do

    do k=2,kbm1
        do j=1,jm
            do i=1,im
                a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.e0*umol)*.5e0 &
                /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
                c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.e0*umol)*.5e0 &
                /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
            end do
        end do
    end do

! the following section solves the equation:
!     dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

! surface and bottom boundary conditions
    const1=(16.6e0**(2.e0/3.e0))*sef

! initialize fields that are not calculated on all boundaries
! but are later used there
    do i=1,im
        do j=1,jm
            l0(i,j)=0.
        end do
    end do
    do i=1,im
        do j=1,jm
            do k=1,kb
                boygr(i,j,k)=0.
                prod(i,j,k)=0.
            end do
        end do
    end do

    do j=1,jmm1
        do i=1,imm1
            utau2(i,j)=sqrt((.5e0*(wusurf(i,j)+wusurf(i+1,j)))**2 &
            +(.5e0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
            uf(i,j,kb)=sqrt((.5e0*(wubot(i,j)+wubot(i+1,j)))**2 &
            +(.5e0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
    end do
    call exchange2d_mpi(utau2,im,jm)
    call exchange2d_mpi(uf(:,:,kb),im,jm)

    do j=1,jm
        do i=1,im
        ! wave breaking energy- a variant of Craig & Banner (1994)
        ! see Mellor and Blumberg, 2003.
            ee(i,j,1)=0.e0
        !          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2(i,j)
            gg(i,j,1)=utau2(i,j)*const1 !KTSIARAS

        ! surface length scale following Stacey (1999).
            l0(i,j)=surfl*utau2(i,j)/grav
        end do
    end do

! calculate speed of sound squared
    do k=1,kb
        do j=1,jm
            do i=1,im
                cc(i,j,k)=0.
            end do
        end do
    end do
    do k=1,kbm1
        do j=1,jm
            do i=1,im
                tp=t(i,j,k)+tbias
                sp=s(i,j,k)+sbias
            ! calculate pressure in units of decibars
                p=grav*rhoref*(-zz(k)*h(i,j))*1.e-4
                cc(i,j,k)=1449.1e0+.00821e0*p+4.55e0*tp-.045e0*tp**2 &
                +1.34e0*(sp-35.0e0)
                cc(i,j,k)=cc(i,j,k) &
                /sqrt((1.e0-.01642e0*p/cc(i,j,k)) &
                *(1.e0-0.40e0*p/cc(i,j,k)**2))
            end do
        end do
    end do

! calculate buoyancy gradient
    do k=2,kbm1
        do j=1,jm
            do i=1,im
                q2b(i,j,k)=abs(q2b(i,j,k))
                q2lb(i,j,k)=abs(q2lb(i,j,k))
                boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k)) &
                /(dzz(k-1)*h(i,j)) &
            ! *** note: comment out next line if dens does not include pressure
                +(grav**2)*2.e0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
            end do
        end do
    end do

    do k=2,kbm1
        do j=1,jm
            do i=1,im
                l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
            !            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
            ! TSIARAS
                gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
                gh(i,j,k)=min(gh(i,j,k),.028e0)
            end do
        end do
    end do

    do j=1,jm
        do i=1,im
        !          l(i,j,1)=kappa*l0(i,j)
            l(i,j,1)=0. !KTSIARAS
            l(i,j,kb)=0.e0
            gh(i,j,1)=0.e0
            gh(i,j,kb)=0.e0
        end do
    end do

! calculate production of turbulent kinetic energy:
    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                prod(i,j,k)=km(i,j,k)*.25e0*sef &
                *((u(i,j,k)-u(i,j,k-1) &
                +u(i+1,j,k)-u(i+1,j,k-1))**2 &
                +(v(i,j,k)-v(i,j,k-1) &
                +v(i,j+1,k)-v(i,j+1,k-1))**2) &
                /(dzz(k-1)*dh(i,j))**2 &
            ! add shear due to internal wave field
                -shiw*km(i,j,k)*boygr(i,j,k)
                prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
            end do
        end do
    end do
    call exchange3d_mpi(prod(:,:,2:kbm1),im,jm,kbm2)

! note: Richardson # dep. dissipation correction (Mellor, 2001; Ezer,
! 2000), depends on ghc the critical number (empirical -6 to -2) to
! increase mixing
    ghc=-6.0e0
    do k=1,kb
        do j=1,jm
            do i=1,im
                stf(i,j,k)=1.e0
            ! It is unclear yet if diss. corr. is needed when surf. waves are included.
                if(gh(i,j,k) < 0.e0) & ! KTSIARAS (commented)
                stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0 !KTSIARAS (commented)
                if(gh(i,j,k) < ghc) stf(i,j,k)=0.1e0 !KTSIARAS (commented)

                dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k) &
                /(b1*l(i,j,k)+small)
            end do
        end do
    end do

    do k=2,kbm1
        do j=1,jm
            do i=1,im
                gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1)) &
                -(2.e0*dti2*dtef(i,j,k)+1.e0))
                ee(i,j,k)=a(i,j,k)*gg(i,j,k)
                gg(i,j,k)=(-2.e0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1) &
                -uf(i,j,k))*gg(i,j,k)
            end do
        end do
    end do

    do k=1,kbm1
        ki=kb-k
        do j=1,jm
            do i=1,im
                uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
            end do
        end do
    end do

! the following section solves the equation:
!     dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
    do j=1,jm
        do i=1,im
            vf(i,j,1)=0.
            vf(i,j,kb)=0.
            ee(i,j,2)=0.e0
            gg(i,j,2)=-kappa*z(2)*dh(i,j)*q2(i,j,2)
            vf(i,j,kb-1)=kappa*(1+z(kbm1))*dh(i,j)*q2(i,j,kbm1)
        end do
    end do
    do k=2,kbm1
        do j=1,jm
            do i=1,im
                dtef(i,j,k)=dtef(i,j,k) &
                *(1.e0+e2*((1.e0/abs(z(k)-z(1)) &
                +1.e0/abs(z(k)-z(kb))) &
                *l(i,j,k)/(dh(i,j)*kappa))**2)
            end do
        end do
    end do
    do k=3,kbm1
        do j=1,jm
            do i=1,im
                gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1)) &
                -(dti2*dtef(i,j,k)+1.e0))
                ee(i,j,k)=a(i,j,k)*gg(i,j,k)
                gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1) &
                +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
            end do
        end do
    end do

    do k=1,kb-2
        ki=kb-k
        do j=1,jm
            do i=1,im
                vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
            end do
        end do
    end do
! the following is to counter the problem of the ratio of two small
! numbers (l = q2l/q2) or one number becoming negative. Two options are
! included below. In this application, the second option, l was less
! noisy when uf or vf is small
    do k=2,kbm1
        do j=1,jm
            do i=1,im
            !           if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
            !             uf(i,j,k)=small
            !             vf(i,j,k)=0.1*dt(i,j)*small
            !           end if
                uf(i,j,k)=abs(uf(i,j,k))
                vf(i,j,k)=abs(vf(i,j,k))
            end do
        end do
    end do

! the following section solves for km and kh
    coef4=18.e0*a1*a1+9.e0*a1*a2
    coef5=9.e0*a1*a2

! note that sm and sh limit to infinity when gh approaches 0.0288
    do k=1,kb
        do j=1,jm
            do i=1,im
                coef1=a2*(1.e0-6.e0*a1/b1*stf(i,j,k))
                coef2=3.e0*a2*b2/stf(i,j,k)+18.e0*a1*a2
                coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf(i,j,k))
                sh(i,j,k)=coef1/(1.e0-coef2*gh(i,j,k))
                sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
                sm(i,j,k)=sm(i,j,k)/(1.e0-coef5*gh(i,j,k))
            end do
        end do
    end do

! there are 2 options for kq which, unlike km and kh, was not derived by
! Mellor and Yamada but was purely empirical based on neutral boundary
! layer data. The choice is whether or not it should be subject to the
! stability factor, sh. Generally, there is not a great difference in
! output
    do k=1,kb
        do j=1,jm
            do i=1,im
                prod(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
                kq(i,j,k)=(prod(i,j,k)*.41e0*sh(i,j,k)+kq(i,j,k))*.5e0
            !            kq(i,j,k)=(prod(i,j,k)*.20+kq(i,j,k))*.5e0
                km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*.5e0
                kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5e0
            end do
        end do
    end do

! cosmetics: make boundr. values as interior (even if not used, printout
! may show strange values)
    do k=1,kb
        do i=1,im
            if(n_north == -1) then
                km(i,jm,k)=km(i,jmm1,k)
                kh(i,jm,k)=kh(i,jmm1,k)
                kq(i,jm,k)=kq(i,jmm1,k)
            end if
            if(n_south == -1) then
                km(i,1,k)=km(i,2,k)
                kh(i,1,k)=kh(i,2,k)
                kq(i,1,k)=kq(i,2,k)
            end if
        end do
        do j=1,jm
            if(n_east == -1) then
                km(im,j,k)=km(imm1,j,k)
                kh(im,j,k)=kh(imm1,j,k)
                kq(im,j,k)=kq(imm1,j,k)
            end if
            if(n_west == -1) then
                km(1,j,k)=km(2,j,k)
                kh(1,j,k)=kh(2,j,k)
                kq(1,j,k)=kq(2,j,k)
            end if
        end do
    end do

    do k=1,kb
        do i=1,im
            do j=1,jm
                km(i,j,k)=km(i,j,k)*fsm(i,j)
                kh(i,j,k)=kh(i,j,k)*fsm(i,j)
                kq(i,j,k)=kq(i,j,k)*fsm(i,j)
            end do
        end do
    end do

    return
    end subroutine profq

! ______________________________________________________________________
    subroutine proft(f,wfsurf,fsurf,nbc)
!!  Solves for vertical diffusion of temperature and salinity using method
!!  described by Richmeyer and Morton (1967)
!!  Note: wfsurf and swrad are negative values when water column is
!!  warming or salt is being added
    implicit none
    include 'pom.h'
    real :: f(im_local,jm_local,kb)
    real :: wfsurf(im_local,jm_local),fsurf(im_local,jm_local)
    integer :: nbc
    real*8 :: a(im,jm,kb),c(im,jm,kb)
    real*8 :: ee(im,jm,kb),gg(im,jm,kb)
    real :: dh(im,jm),rad(im,jm,kb)
    real :: r(5),ad1(5),ad2(5)
    integer :: i,j,k,ki

! irradiance parameters after Paulson and Simpson (1977)
!       ntp               1      2       3       4       5
!   Jerlov type           i      ia      ib      ii     iii
    data r   /       .58e0,  .62e0,  .67e0,  .77e0,  .78e0 /
    data ad1 /       .35e0,  .60e0,  1.0e0,  1.5e0,  1.4e0 /
    data ad2 /       23.e0,  20.e0,  17.e0,  14.e0,  7.9e0 /

! surface boundary condition:
!       nbc   prescribed    prescribed   short wave
!             temperature      flux      penetration
!             or salinity               (temperature
!                                           only)
!        1        no           yes           no
!        2        no           yes           yes
!        3        yes          no            no
!        4        yes          no            yes
! note that only 1 and 3 are allowed for salinity

! the following section solves the equation
!     dti2*(kh*f')'-f=-fb

    do j=1,jm
        do i=1,im
            dh(i,j)=h(i,j)+etf(i,j)
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                a(i,j,k)=0.e0
                c(i,j,k)=0.e0
                ee(i,j,k)=0.e0
                gg(i,j,k)=0.e0
            end do
        end do
    end do

    do k=2,kbm1
        do j=1,jm
            do i=1,im
                a(i,j,k-1)=-dti2*(kh(i,j,k)+umol) &
                /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
                c(i,j,k)=-dti2*(kh(i,j,k)+umol) &
                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
            end do
        end do
    end do

! calculate penetrative radiation. At the bottom any unattenuated
! radiation is deposited in the bottom layer
    do k=1,kb
        do j=1,jm
            do i=1,im
                rad(i,j,k)=0.e0
            end do
        end do
    end do

    if(nbc == 2 .OR. nbc == 4) then
        do k=1,kbm1
            do j=1,jm
                do i=1,im
                    rad(i,j,k)=swrad(i,j) &
                    *(r(ntp)*exp(z(k)*dh(i,j)/ad1(ntp)) &
                    +(1.e0-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)))
                end do
            end do
        end do
    end if

    if(nbc == 1) then

        do j=1,jm
            do i=1,im
                ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
                gg(i,j,1)=dti2*wfsurf(i,j)/(dz(1)*dh(i,j))-f(i,j,1)
                gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
            end do
        end do

    else if(nbc == 2) then

        do j=1,jm
            do i=1,im
                ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
                gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2)) &
                /(dz(1)*dh(i,j)) &
                -f(i,j,1)
                gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
            end do
        end do

    else if(nbc == 3 .OR. nbc == 4) then

        do j=1,jm
            do i=1,im
                ee(i,j,1)=0.e0
                gg(i,j,1)=fsurf(i,j)
            end do
        end do

    end if

    do k=2,kbm2
        do j=1,jm
            do i=1,im
                gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
                ee(i,j,k)=a(i,j,k)*gg(i,j,k)
                gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k) &
                +dti2*(rad(i,j,k)-rad(i,j,k+1)) &
                /(dh(i,j)*dz(k))) &
                *gg(i,j,k)
            end do
        end do
    end do

! bottom adiabatic boundary condition
    do j=1,jm
        do i=1,im
            if(fsm(i,j) == 1.)then
                f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1) &
                +dti2*(rad(i,j,kbm1)-rad(i,j,kb)) &
                /(dh(i,j)*dz(kbm1))) &
                /(c(i,j,kbm1)*(1.e0-ee(i,j,kbm2))-1.e0)
            endif
        end do
    end do
     

    do k=2,kbm1
        ki=kb-k
        do j=1,jm
            do i=1,im
                if(fsm(i,j) == 1.)then
                    f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
                endif
            end do
        end do
    end do

    return
    end subroutine proft

! ______________________________________________________________________
    subroutine profu
! Solves for vertical diffusion of x-momentum using method described by
! Richmeyer and Morton (1967)
! Note: wusurf has the opposite sign to the wind speed
    implicit none
    include 'pom.h'
    real*8 :: a(im,jm,kb),c(im,jm,kb)
    real*8 :: ee(im,jm,kb),gg(im,jm,kb)
    real :: dh(im,jm)
    integer :: i,j,k,ki

! the following section solves the equation
!   dti2*(km*u')'-u=-ub
    do j=1,jm
        do i=1,im
            dh(i,j)=1.e0
        end do
    end do

    do j=2,jm
        do i=2,im
            dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5e0
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                a(i,j,k)=0.e0
                c(i,j,k)=0.e0
                ee(i,j,k)=0.e0
                gg(i,j,k)=0.e0
            end do
        end do
    end do

    do k=1,kb
        do j=2,jm
            do i=2,im
                c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5e0
            end do
        end do
    end do

    do k=2,kbm1
        do j=1,jm
            do i=1,im
                a(i,j,k-1)=-dti2*(c(i,j,k)+umol) &
                /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
                c(i,j,k)=-dti2*(c(i,j,k)+umol) &
                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
            end do
        end do
    end do

    do j=1,jm
        do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j)) &
            -uf(i,j,1)) &
            /(a(i,j,1)-1.e0)
        end do
    end do

    do k=2,kbm2
        do j=1,jm
            do i=1,im
                gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
                ee(i,j,k)=a(i,j,k)*gg(i,j,k)
                gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
            end do
        end do
    end do

    do j=2,jmm1
        do i=2,imm1
            if(dum(i,j) == 1.)then !KTSIARAS to remove NaNs in land points
                tps(i,j)=0.5e0*(cbc(i,j)+cbc(i-1,j)) &
                *sqrt(ub(i,j,kbm1)**2 &
                +(.25e0*(vb(i,j,kbm1)+vb(i,j+1,kbm1) &
                +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
                uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1)) &
                /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0 &
                -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
                uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
            endif
        end do
    end do

    do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
            do i=2,imm1
                if(dum(i,j) == 1.)then !KTSIARAS to remove NaNs in land points
                    uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
                endif
            end do
        end do
    end do

    do j=2,jmm1
        do i=2,imm1
            wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
    end do
    call exchange2d_mpi(wubot,im,jm)

    return
    end subroutine profu

! ______________________________________________________________________
    subroutine profv
!!  Solves for vertical diffusion of x-momentum using method described by
!!  Richmeyer and Morton (1967)
!!  Note: wvsurf has the opposite sign to the wind speed
    implicit none
    include 'pom.h'
    real*8 :: a(im,jm,kb),c(im,jm,kb)
    real*8 :: ee(im,jm,kb),gg(im,jm,kb)
    real :: dh(im,jm)
    integer :: i,j,k,ki

! the following section solves the equation
!     dti2*(km*u')'-u=-ub

    do j=1,jm
        do i=1,im
            dh(i,j)=1.e0
        end do
    end do

    do j=2,jm
        do i=2,im
            dh(i,j)=.5e0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
    end do

    do k=1,kb
        do j=1,jm
            do i=1,im
                a(i,j,k)=0.e0
                c(i,j,k)=0.e0
                ee(i,j,k)=0.e0
                gg(i,j,k)=0.e0
            end do
        end do
    end do

    do k=1,kb
        do j=2,jm
            do i=2,im
                c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5e0
            end do
        end do
    end do

    do k=2,kbm1
        do j=1,jm
            do i=1,im
                a(i,j,k-1)=-dti2*(c(i,j,k)+umol) &
                /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
                c(i,j,k)=-dti2*(c(i,j,k)+umol) &
                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
            end do
        end do
    end do

    do j=1,jm
        do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1)) &
            /(a(i,j,1)-1.e0)
        end do
    end do

    do k=2,kbm2
        do j=1,jm
            do i=1,im
                gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
                ee(i,j,k)=a(i,j,k)*gg(i,j,k)
                gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
            end do
        end do
    end do

    do j=2,jmm1
        do i=2,imm1
            if(dvm(i,j) == 1.)then !KTSIARAS to remove NaNs in land points
                tps(i,j)=0.5e0*(cbc(i,j)+cbc(i,j-1)) &
                *sqrt((.25e0*(ub(i,j,kbm1)+ub(i+1,j,kbm1) &
                +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2 &
                +vb(i,j,kbm1)**2)
                vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1)) &
                /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0 &
                -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
                vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
            endif
        end do
    end do

    do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
            do i=2,imm1
                if(dvm(i,j) == 1.)then !KTSIARAS to remove NaNs in land points
                    vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
                endif
            end do
        end do
    end do

    do j=2,jmm1
        do i=2,imm1
            wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
    end do
    call exchange2d_mpi(wvbot,im,jm)

    return
    end subroutine profv

! ______________________________________________________________________
    subroutine smol_adif(xmassflux,ymassflux,zwflux,ff)
!!  Calculate the antidiffusive velocity used to reduce the numerical
!!  diffusion associated with the upstream differencing scheme
!!  this is based on a subroutine of Gianmaria Sannino (Inter-university
!!  Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
!!  National Agency for New Technology and Environment, Rome, Italy)
    implicit none
    include 'pom.h'
    real :: ff(im_local,jm_local,kb)
    real :: xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
    real :: mol,abs_1,abs_2
    real :: value_min,epsilon
    real :: udx,u2dt,vdy,v2dt,wdz,w2dt
    integer :: i,j,k
    parameter (value_min=1.e-9,epsilon=1.0e-14)

! apply temperature and salinity mask
    do k=1,kb
        do i=1,im
            do j=1,jm
                ff(i,j,k)=ff(i,j,k)*fsm(i,j)
            end do
        end do
    end do

! recalculate mass fluxes with antidiffusion velocity
    do k=1,kbm1
        do j=2,jmm1
            do i=2,im
                if(ff(i,j,k) < value_min .OR. &
                ff(i-1,j,k) < value_min) then
                    xmassflux(i,j,k)=0.e0
                else
                    udx=abs(xmassflux(i,j,k))
                    u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.e0 &
                    /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
                    mol=(ff(i,j,k)-ff(i-1,j,k)) &
                    /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
                    xmassflux(i,j,k)=(udx-u2dt)*mol*sw
                    abs_1=abs(udx)
                    abs_2=abs(u2dt)
                    if(abs_1 < abs_2) xmassflux(i,j,k)=0.e0
                end if
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jm
            do i=2,imm1
                if(ff(i,j,k) < value_min .OR. &
                ff(i,j-1,k) < value_min) then
                    ymassflux(i,j,k)=0.e0
                else
                    vdy=abs(ymassflux(i,j,k))
                    v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.e0 &
                    /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
                    mol=(ff(i,j,k)-ff(i,j-1,k)) &
                    /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
                    ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
                    abs_1=abs(vdy)
                    abs_2=abs(v2dt)
                    if(abs_1 < abs_2) ymassflux(i,j,k)=0.e0
                end if
            end do
        end do
    end do

    do k=2,kbm1
        do j=2,jmm1
            do i=2,imm1
                if(ff(i,j,k) < value_min .OR. &
                ff(i,j,k-1) < value_min) then
                    zwflux(i,j,k)=0.e0
                else
                    wdz=abs(zwflux(i,j,k))
                    w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/ &
                    (dzz(k-1)*dt(i,j))
                    mol=(ff(i,j,k-1)-ff(i,j,k)) &
                    /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
                    zwflux(i,j,k)=(wdz-w2dt)*mol*sw
                    abs_1=abs(wdz)
                    abs_2=abs(w2dt)
                    if(abs_1 < abs_2)zwflux(i,j,k)=0.e0
                end if
            end do
        end do
    end do

    return
    end subroutine smol_adif

! ______________________________________________________________________
    subroutine vertvl
!!  Calculates vertical velocity
    implicit none
    include 'pom.h'
    real :: xflux(im,jm,kb),yflux(im,jm,kb)
    integer :: i,j,k

    do k=1,kb
        do j=1,jm
            do i=1,im
                xflux(i,j,k)=0.e0
                yflux(i,j,k)=0.e0
            end do
        end do
    end do

! reestablish boundary conditions
    do k=1,kbm1
        do j=2,jm
            do i=2,im
                xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)) &
                *(dt(i,j)+dt(i-1,j))*u(i,j,k)
            end do
        end do
    end do

    do k=1,kbm1
        do j=2,jm
            do i=2,im
                yflux(i,j,k)=.25e0*(dx(i,j)+dx(i,j-1)) &
                *(dt(i,j)+dt(i,j-1))*v(i,j,k)
            end do
        end do
    end do

! note: if one wishes to include freshwater flux, the surface velocity
! should be set to vflux(i,j). See also change made to 2-D volume
! conservation equation which calculates elf
    do j=2,jmm1
        do i=2,imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
        end do
    end do

    do k=1,kbm1
        do j=2,jmm1
            do i=2,imm1
                w(i,j,k+1)=w(i,j,k) &
                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k) &
                +yflux(i,j+1,k)-yflux(i,j,k)) &
                /(dx(i,j)*dy(i,j)) &
                +(etf(i,j)-etb(i,j))/dti2)
            end do
        end do
    end do

    return
    end subroutine vertvl

! ______________________________________________________________________
subroutine realvertvl
!!  Calculates real vertical velocity (wr)
    implicit none
    include 'pom.h'
    integer :: i,j,k
    real :: dxr,dxl,dyt,dyb
    do k=1,kb
        do j=1,jm
            do i=1,im
                wr(i,j,k)=0.
            end do
        end do
    end do

    do k=1,kbm1
        do j=1,jm
            do i=1,im
                tps(i,j)=zz(k)*dt(i,j) + et(i,j)
            end do
        end do
        do j=2,jmm1
            do i=2,imm1
                dxr=2.0/(dx(i+1,j)+dx(i,j))
                dxl=2.0/(dx(i,j)+dx(i-1,j))
                dyt=2.0/(dy(i,j+1)+dy(i,j))
                dyb=2.0/(dy(i,j)+dy(i,j-1))
                wr(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))+0.5* &
                (u(i+1,j,k)*(tps(i+1,j)-tps(i,j))*dxr+ &
                u(i,j,k)*(tps(i,j)-tps(i-1,j))*dxl+ &
                v(i,j+1,k)*(tps(i,j+1)-tps(i,j))*dyt+ &
                v(i,j,k)*(tps(i,j)-tps(i,j-1))*dyb) &
                +(1.0+zz(k))*(etf(i,j)-etb(i,j))/dti2
            end do
        end do
    end do

    call exchange3d_mpi(wr(:,:,1:kbm1),im,jm,kbm1)

    do k=1,kb
        do i=1,im
            if(n_south == -1) wr(i,1,k)=wr(i,2,k)
            if(n_north == -1) wr(i,jm,k)=wr(i,jmm1,k)
        end do
    end do
    do k=1,kb
        do j=1,jm
            if(n_west == -1) wr(1,j,k)=wr(2,j,k)
            if(n_east == -1) wr(im,j,k)=wr(imm1,j,k)
        end do
    end do

    do k=1,kbm1
        do j=1,jm
            do i=1,im
                wr(i,j,k)=fsm(i,j)*wr(i,j,k)
            end do
        end do
    end do

    return
    end subroutine realvertvl




