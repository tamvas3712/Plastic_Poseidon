# 1 "more_subs.F90"
integer function randomf(n1,n2)
    implicit none
    integer :: n1,n2
    real*4 :: r
    call random_number(r)
    randomf=n1+int((n2-n1+1)*r)
end function randomf


subroutine init_lagr
    include 'pom.h'
    include 'lagr.h'

!    calculate initial cartesian coordinates (xo,yo) from Lon/Lat (xsi,ysi)
!    and set initial depth (fdep) in the water column
    do k=1,nsi
        if(variant(k) /= 0)then
            call Geo2Cart (xo(k), yo(k), alon(1,1)-(stepx/2.), &
            alat(1,1)-(stepy/2.), xsi(k),ysi(k), pi, dg)
            call FCoord (io, jo, alon(1,1)-(stepx/2.), &
            alat(1,1)-(stepy/2.), stepx, stepy, &
            xo(k), yo(k), pi, dg)
            zo(k)=fdep(k)
!        zo(k)=-10.
        endif
    enddo

!   set x-max and y-max (used to check if inside model domain)
    xm = 0.
    do i = 1, im
        xm = xm + dx(i,1)
    end do
    ym = 0.
    do j = 1, jm
        ym = ym + dy(1,j)
    end do

! Initialize lbc (find different boundary conditions on land boxes)
    call land_bc(lbc, h, im, jm, 1.)

! Assign retention time for coastal areas
    do j = 1, jm
        do i = 1, im
            rbc(i, j) = 0.
!            if(h(i, j) .le. 1.) rbc(i, j) = rtime(nsp,ivariant)
            if(h(i, j) <= 1.) rbc(i, j) = 0.
        end do
    end do

    return
end subroutine init_lagr

subroutine totfish(iwrite)
!! calculates total #SIs, biomass/particles & mean depth/density/weight
    include 'ibmall.fcm'
    integer :: ch,i,nn,iwres
    integer :: icountempty,iwrite
    real :: weira(nvariantsmax,nspecies)

    do nsp=1,nspecies
        do i=1,nvariants(nsp)
            weivariant(i,nsp)=0.
            denvariant(i,nsp)=0.
            depvariant(i,nsp)=0.
            tbvariant(i,nsp)=0.
            tnvariant(i,nsp)=0.
            tsivariant(i,nsp)=0.
            tbeach(i,nsp)=0.
        enddo
    enddo

! Calculation of empty SIs,then fill with new eggs in ersem_fishGrad.F
    do n=1,nmax
        nempty(n)=0
    enddo
    icountempty=0

    do n=1,nsi
        if(variant(n) == 0)then
            icountempty=icountempty+1
            nempty(icountempty)=n
        else
            nsp=species(n)
            ivariant=variant(n)

            tbvariant(ivariant,nsp)=tbvariant(ivariant,nsp)+nfi(n)*xfi(n)
            tnvariant(ivariant,nsp)=tnvariant(ivariant,nsp)+nfi(n)
            tsivariant(ivariant,nsp)=tsivariant(ivariant,nsp)+1
            if(ref(n) > 0.)tbeach(ivariant,nsp)=tbeach(ivariant,nsp)+1
            weivariant(ivariant,nsp)=weivariant(ivariant,nsp)+nfi(n)*xfi(n)
            depvariant(ivariant,nsp)=depvariant(ivariant,nsp)+nfi(n)*fdep(n)
        endif
    enddo

    do n=1,nspecies
        do i=1,nvariants(n)
            if(tnvariant(i,n) > 0.)then
                weivariant(i,n)=weivariant(i,n)/tnvariant(i,n)
                depvariant(i,n)=depvariant(i,n)/tnvariant(i,n)
            else
                weivariant(i,n)=0.
                depvariant(i,n)=0.
            endif
        enddo
    enddo

    if(iwrite == 1)then
        write(*,*)'NSI=',nsi
        write(*,*)'EMPTY=',icountempty
        do n=1,nspecies
            write(*,*)'SPECIES',n
            adultbiom=0.
            do ivariant=1,nvariants(n)
                write(*,*)'TOTAL BIOMASS',ivariant,tbvariant(ivariant,n)*1.e-6
                write(*,*)'PARTICLES',ivariant,tnvariant(ivariant,n)
                write(*,*)'DENSITY',ivariant,denvariant(ivariant,n)
                write(*,*)'DEPTH',ivariant,depvariant(ivariant,n)
                write(*,*)'SIs',ivariant,tsivariant(ivariant,n)
                write(*,*)'SIs-BEACHED',ivariant,tbeach(ivariant,n)
                write(*,*)'PARTICLE WEIGHT(mgr)',ivariant,weivariant(ivariant,n)
            enddo
        enddo
    endif

    return
end subroutine totfish



