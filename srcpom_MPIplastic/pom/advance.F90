subroutine advance
!!  advance POM 1 step in time
    implicit none
    include 'pom.h'
    real :: delt,dti_fish
    integer :: n,i
! get time
    call get_time

! update hydrodynamics
    call update_hydro
    write(*,*)'update_hydro'
    #ifdef drift
    call windwave
    write(*,*)'windwave'
    #endif

   call lagr_move(dti) !lagrangian model 

    if(mod(iint,irestart) == 0)then
        call write_ibmrestart_pnetcdf !write ibm output (netcdf)
        call totfish(1)               !write diagnostics on output     
        write(*,*)'writeibm'
    endif

    return
end subroutine advance

