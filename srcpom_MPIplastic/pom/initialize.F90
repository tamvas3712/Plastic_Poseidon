    subroutine initialize
!!  initialize Model: 
!!  1) reads grid, hydro fields, initial IBM and date/timestep
!!  2) initializes MPI(parallel)

    implicit none
    include 'pom.h'
    integer :: i,j,k,iyr_n,ihour
    real ::    year               !KTSIARAS

! initialize the MPI execution environment and create communicator for
! Internal POM communications
    call initialize_mpi

! distribute the model domain across processors
    call distribute_mpi_ibm

! read time step, write frequency and some constants
    call read_input

! read model grid
    call read_grid

    write(*,*)'GRID',alon(1,1),alat(1,1)
!--read time from current date -------------------------------
    open(15,file='dir1')
    read(15,'(i2.2,i2.2,i2.2)') iday,imon,iyr1
    close(15)

    if(iyr1 <= 50)then
        iyr=iyr1+2000
    else
        iyr=iyr1+1900
    endif

    year=float(iyr)
    IHOUR = 18
    call findtime0(time0,year,ihour,iday,imon,iyr)
    TIME = TIME0
    write(*,*)'TEST-TIME',time0,iday,imon,iyr

! read hydrodynamics
    call read_hydroibm_pnetcdf
    write(*,*)'OK'
! update hydrodynamics
    call update_hydro
    call density(s,t,rho)

    write(*,*)'OK'
!-read plastics
    call read_ibm_pnetcdf
!- initialize lagrangian
    call init_lagr
!--read bacteria 3-d field (for biofouling) / from ERSEM simulation
    call read_bac3d_pnetcdf

!------------------------------------------------------------
! check for errors
    call sum0d_mpi(error_status,master_task)
    call bcast0d_mpi(error_status,master_task)
    if(error_status /= 0) then
        if(my_task == master_task) write(*,'(/a)') &
        'POM terminated with error'
        call finalize_mpi
        stop
    end if

    if(my_task == master_task) write(*,'(/a)') 'End of initialization'

    return
    end subroutine initialize


subroutine read_input
!!  Read input values and define constants
    implicit none
    include 'pom.h'
    namelist/pom_nml/ title,netcdf_file,dti, &
    time_start,write_rst, &
    write_rst_file,days,prtd1,prtd2, &
    swtch

! read input namelist (time step, time to write etc)
    open(73,file='pomibm.nml',status='old')
    read(73,nml=pom_nml)
    close(73)

! Input of filenames and constants

! gravity constant (S.I. units)
    grav=9.806e0

! End of input of constants

! calculate some constants
    small=1.e-9           ! Small value
    pi=atan(1.e0)*4.e0    ! PI

    dti2=dti*2

    iend=max0(nint(days*24.e0*3600.e0/dti),2)
    iprint=nint(prtd1*24.e0*3600.e0/dti)
    iswtch=nint(swtch*24.e0*3600.e0/dti)
    irestart=nint(write_rst*24.e0*3600.e0/dti)
     

! initialise time
    time0=0.e0
    time=0.e0

! print initial summary
    if(my_task == master_task) then
        write(6,'(/'' title      = '',a40)') title
        write(6,'('' dti        = '',f10.1)') dti
        write(6,'('' time_start = '',a26)') time_start
        write(6,'('' days       = '',f10.4)') days
        write(6,'('' iend       = '',i10)') iend
        write(6,'('' prtd1      = '',f10.4)') prtd1
        write(6,'('' iprint     = '',i10)') iprint
        write(6,'('' prtd2      = '',f10.4)') prtd2
        write(6,'('' swtch      = '',f10.2)') swtch
        write(6,'('' iswtch     = '',i10)') iswtch
        write(6,'('' grav       = '',f10.4)') grav
    end if

    return
    end subroutine read_input

subroutine read_grid
!!  Set up vertical and horizontal grid, topography, areas and masks
    implicit none
    include 'pom.h'
    integer :: i,j,k
    real :: deg2rad

! degrees to radians
    deg2rad=pi/180.

! read grid
    call read_grid_pnetcdf

    do j=1,jm
        do i=1,im
            alon(i,j)=east_e(i,j)
            alat(i,j)=north_e(i,j)
        enddo
    enddo

! derived vertical grid variables
    do k=1,kb-1
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
    end do
    dz(kb)=0.
    dzz(kb)=0.

! print vertical grid information
    if(my_task == master_task) then
        write(6,'(/2x,a,7x,a,9x,a,9x,a,9x,a)') 'k','z','zz','dz','dzz'
        do k=1,kb
            write(6,'(1x,i5,4f10.3)') k,z(k),zz(k),dz(k),dzz(k)
        end do
    end if

    return
    end subroutine read_grid


subroutine update_hydro ! KTSIARAS / for Off-line Coupling
!! calculates at every time-step the hydro fields with linear interpolation between 3-hour fields
    implicit none
    include 'pom.h'
    include 'hydro.h'

    integer :: nstart
    parameter (nstart=0)

    integer :: i,j,k,is
    real :: xp,xm,factor,tday,thour
    integer :: itday,mp,mm,ntime
    real :: xflux(im,jm,kb),yflux(im,jm,kb)
    thour=nfre ! time between wind files (hours)
    tday=(nfre*3600./86400.) !(days)
    itday=int(thour*3600.e0/dti)

! linear interpolation in time
    ntime=(time-time0)/tday
    mp=ntime-nstart+1
    mm=mp+1
    xm=(time-time0)/tday-ntime
    xp=1.-xm

!    write(*,*)'TESTCLOCK',time,mp,mm,xp,xm

! linear interpolation in time

    do j=1,jm
        do i=1,im
            elb(i,j)=(el_f(i,j,mp)*xp+el_f(i,j,mm)*xm)*fsm(i,j)
            dt(i,j)=elb(i,j)+h(i,j)
            do k=1,kb-1
                u(i,j,k)=(u_f(i,j,k,mp)*xp+u_f(i,j,k,mm)*xm)*dum(i,j)
                v(i,j,k)=(v_f(i,j,k,mp)*xp+v_f(i,j,k,mm)*xm)*dvm(i,j)
                w(i,j,k)=(w_f(i,j,k,mp)*xp+w_f(i,j,k,mm)*xm)*fsm(i,j)
                t(i,j,k)=(t_f(i,j,k,mp)*xp+t_f(i,j,k,mm)*xm)*fsm(i,j)
                s(i,j,k)=(s_f(i,j,k,mp)*xp+s_f(i,j,k,mm)*xm)*fsm(i,j)
                kh(i,j,k)=(kh_f(i,j,k,mp)*xp+kh_f(i,j,k,mm)*xm)*fsm(i,j)
                aam(i,j,k)=(aam_f(i,j,k,mp)*xp+aam_f(i,j,k,mm)*xm)* &
                fsm(i,j)
            end do
        enddo
    enddo

    return
    end subroutine update_hydro


    subroutine windwave !KTSIARAS
!!  Read and interpolate (in time) wind and wave
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

! read initial wind/wave files
    if (iint == 1) then
        call read_wave_pnetcdf(iint/iflux+1) 
        call read_wind2_pnetcdf(iint/iflux+1)
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
            uair(i,j)=fold*uairf(i,j,iold)+fnew*uairf(i,j,inew)           !U-wind at 10m above surface
            vair(i,j)=fold*vairf(i,j,iold)+fnew*vairf(i,j,inew)           !V-wind at 10m above surface
            wavehei(i,j)=fold*waveheif(i,j,iold)+fnew*waveheif(i,j,inew)  !Significant wave height
            waveper(i,j)=fold*waveperf(i,j,iold)+fnew*waveperf(i,j,inew)  !wave period
            stokesx(i,j)=fold*stokesxf(i,j,iold)+fnew*stokesxf(i,j,inew)  !stokes drift-x
            stokesy(i,j)=fold*stokesyf(i,j,iold)+fnew*stokesyf(i,j,inew)  !stokes drift-y
        end do
    end do

    write(*,*)'TEST-INTERP',time,my_task,ntime,tfluxd,fnew,fold,inew, &
    iold,waveper(329,150),wavehei(329,150),stokesx(329,150)

    return
    end subroutine windwave

