! parallel_mpi.f

! subroutines for communicating between processors using MPI

! ______________________________________________________________________
    subroutine initialize_mpi
!!  Set up MPI execution environment and define the POM communicator
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: ierr
! initiate MPI environment
    call mpi_init(ierr)
! determine processor rank
    call mpi_comm_rank(mpi_comm_world,my_task,ierr)
    pom_comm=mpi_comm_world
    master_task=0
    error_status=0
    return
    end subroutine initialize_mpi

! ______________________________________________________________________
    subroutine finalize_mpi
!!  Terminate the MPI execution environment
    implicit none
    include 'mpif.h'
    integer :: ierr
! terminate MPI environment
    call mpi_finalize(ierr)
    return
    end subroutine finalize_mpi

! ______________________________________________________________________
    subroutine distribute_mpi
!!  Distribute the model domain across processors
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: i,j,ierr,nproc,nproc_x,nproc_y

! determine the number of processors
    call mpi_comm_size(pom_comm,nproc,ierr)

! check number of processors
    if(nproc /= n_proc) then
        error_status=1
        if(my_task == master_task) write(*,'(a//a)') &
        'Incompatible number of processors','POM terminated with error'
        call finalize_mpi
        stop
    end if

! determine the number of processors in x
    if(mod(im_global-2,im_local-2) == 0) then
        nproc_x=(im_global-2)/(im_local-2)
    else
        nproc_x=(im_global-2)/(im_local-2) + 1
    end if

! determine the number of processors in y
    if(mod(jm_global-2,jm_local-2) == 0) then
        nproc_y=(jm_global-2)/(jm_local-2)
    else
        nproc_y=(jm_global-2)/(jm_local-2) + 1
    end if

    write(*,*)'CPU',nproc_x,nproc_y,nproc_x*nproc_y,n_proc
! check local size
    if(nproc_x*nproc_y > n_proc) then
        error_status=1
        if(my_task == master_task) write(*,'(a//a)') &
        'im_local or jm_local is too low','POM terminated with error'
        call finalize_mpi
        stop
    end if

! detemine global and local indices
    im=im_local
    do i=1,im_local
        i_global(i)=0
    end do
    do i=1,im
        i_global(i)=i+mod(my_task,nproc_x)*(im-2)
        if(i_global(i) > im_global) then
            im=i-1
            i_global(i)=0
            cycle
        end if
    end do
    imm1=im-1
    imm2=im-2

    jm=jm_local
    do j=1,jm_local
        j_global(j)=0
    end do
    do j=1,jm
        j_global(j)=j+(my_task/nproc_x)*(jm-2)
        if(j_global(j) > jm_global) then
            jm=j-1
            j_global(j)=0
            cycle
        end if
    end do
    jmm1=jm-1
    jmm2=jm-2

    kbm1=kb-1
    kbm2=kb-2

! determine the neighbors (tasks)
    n_east=my_task+1
    n_west=my_task-1
    n_north=my_task+nproc_x
    n_south=my_task-nproc_x

    if(mod(n_east,nproc_x) == 0) n_east=-1
    if(mod(n_west+1,nproc_x) == 0) n_west=-1
    if(n_north/nproc_x == nproc_y) n_north=-1
    if((n_south+nproc_x)/nproc_x == 0) n_south=-1

    write(*,*)'TEST INDICES',my_task,i_global(1),i_global(im), &
    j_global(1),j_global(jm), &
    n_east,n_west,n_north,n_south

    return
    end subroutine distribute_mpi

    


subroutine distribute_mpi_ibm
!!  Distribute the model domain across processors
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: i,j,ierr,nproc,nproc_x,nproc_y
    character(120) :: netcdf_ic_file
    integer :: ncid,status,nsi_id
    integer :: vdims(4)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start(1),edge(1)

! determine the number of processors
    call mpi_comm_size(pom_comm,nproc,ierr)

! check number of processors
    if(nproc /= n_proc) then
        error_status=1
        if(my_task == master_task) write(*,'(a//a)') &
        'Incompatible number of processors','POM terminated with error'
        call finalize_mpi
        stop
    end if

    nproc_x=1
    nproc_y=1

    n_east=-1
    n_west=-1
    n_north=-1
    n_south=-1

    im=im_global
    do i=1,im
        i_global(i)=i
    end do
    imm1=im-1
    imm2=im-2

    jm=jm_global
    do j=1,jm
        j_global(j)=j
    end do
    jmm1=jm-1
    jmm2=jm-2


    kbm1=kb-1
    kbm2=kb-2

! detemine global and local indices (should be nsi_global/nproc=integer)

! check number of processors
    if(mod(nsi_global,nproc) /= 0) then
        error_status=1
        if(my_task == master_task) write(*,'(a//a)') &
        'ERROR NSI/number of processors should be integer '
        call finalize_mpi
        stop
    end if

    nsi=nsi_local

    do i=1,nsi
        n_global(i)=0
    end do
    do i=1,nsi
        n_global(i)=i+mod(my_task,nproc)*nsi
        if(n_global(i) > nsi_global) then
            nsi=i-1
            n_global(i)=0
            cycle
        end if
    end do

!      nmaxeggs=nsiVariants(1,1)
!      do i=1,nmaxeggs
!       negg_global(i)=nsi_global+(i-1)*nproc+1+mod(my_task,nproc)
!        negg_global(i)=nsi_globalNew+i+mod(my_task,nproc)*nmaxeggs
!      enddo

    write(*,*)'TEST INDICES',my_task,nsi,n_global(1),n_global(nsi)

    return
    end subroutine distribute_mpi_ibm


    

subroutine sum0d_mpi(work,to)
!!  Send sum of WORK to node TO
    implicit none
    integer :: to
    real :: work,tmp
    include 'mpif.h'
    include 'pom.h'
    integer :: ierr
! sum data
    call mpi_reduce(work,tmp,1,mpi_real,mpi_sum,to,pom_comm,ierr)
    work=tmp
    return
    end subroutine sum0d_mpi

    

subroutine bcast0d_mpi(work,from)
!!  Send WORK to all nodes from node FROM
    implicit none
    integer :: from
    real :: work
    include 'mpif.h'
    include 'pom.h'
    integer :: ierr
! broadcast data
    call mpi_bcast(work,1,mpi_real,from,pom_comm,ierr)
    return
    end subroutine bcast0d_mpi

    


subroutine exchange2d_mpi(work,nx,ny)
!!  Exchange ghost cells around 2d local grids
!!  one band at a time
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: nx,ny
    real :: work(nx,ny)
    integer :: i,j,k
    integer :: ierr
    integer :: istatus(mpi_status_size)
    real :: send_east(ny),recv_west(ny)
    real :: send_west(ny),recv_east(ny)
    real :: send_north(nx),recv_south(nx)
    real :: send_south(nx),recv_north(nx)

! send ghost cell data to the east
    if(n_east /= -1) then
        do j=1,ny
            send_east(j)=work(nx-1,j)
        end do
        call mpi_send(send_east,ny,mpi_real,n_east,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the west
    if(n_west /= -1) then
        call mpi_recv(recv_west,ny,mpi_real,n_west,n_west, &
        pom_comm,istatus,ierr)
        do j=1,ny
            work(1,j)=recv_west(j)
        end do
    end if

! send ghost cell data to the west
    if(n_west /= -1) then
        do j=1,ny
            send_west(j)=work(2,j)
        end do
        call mpi_send(send_west,ny,mpi_real,n_west,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the east
    if(n_east /= -1) then
        call mpi_recv(recv_east,ny,mpi_real,n_east,n_east, &
        pom_comm,istatus,ierr)
        do j=1,ny
            work(nx,j)=recv_east(j)
        end do
    end if

! send ghost cell data to the north
    if(n_north /= -1) then
        do i=1,nx
            send_north(i)=work(i,ny-1)
        end do
        call mpi_send(send_north,nx,mpi_real,n_north,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the south
    if(n_south /= -1) then
        call mpi_recv(recv_south,nx,mpi_real,n_south,n_south, &
        pom_comm,istatus,ierr)
        do i=1,nx
            work(i,1)=recv_south(i)
        end do
    end if

! send ghost cell data to the south
    if(n_south /= -1) then
        do i=1,nx
            send_south(i)=work(i,2)
        end do
        call mpi_send(send_south,nx,mpi_real,n_south,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the north
    if(n_north /= -1) then
        call mpi_recv(recv_north,nx,mpi_real,n_north,n_north, &
        pom_comm,istatus,ierr)
        do i=1,nx
            work(i,ny)=recv_north(i)
        end do
    end if

    return
    end subroutine exchange2d_mpi

    

subroutine exchange3d_mpi(work,nx,ny,nz)
!!  Exchange ghost cells around 3d local grids
!!  one band at a time
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: nx,ny,nz
    real :: work(nx,ny,nz)
    integer :: i,j,k
    integer :: ierr
    integer :: istatus(mpi_status_size)
    real :: send_east(ny*nz),recv_west(ny*nz)
    real :: send_west(ny*nz),recv_east(ny*nz)
    real :: send_north(nx*nz),recv_south(nx*nz)
    real :: send_south(nx*nz),recv_north(nx*nz)

! send ghost cell data to the east
    if(n_east /= -1) then
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                send_east(i)=work(nx-1,j,k)
            end do
        end do
        call mpi_send(send_east,ny*nz,mpi_real,n_east,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the west
    if(n_west /= -1) then
        call mpi_recv(recv_west,ny*nz,mpi_real,n_west,n_west, &
        pom_comm,istatus,ierr)
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                work(1,j,k)=recv_west(i)
            end do
        end do
    end if

! send ghost cell data to the west
    if(n_west /= -1) then
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                send_west(i)=work(2,j,k)
            end do
        end do
        call mpi_send(send_west,ny*nz,mpi_real,n_west,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the east
    if(n_east /= -1) then
        call mpi_recv(recv_east,ny*nz,mpi_real,n_east,n_east, &
        pom_comm,istatus,ierr)
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                work(nx,j,k)=recv_east(i)
            end do
        end do
    end if

! send ghost cell data to the north
    if(n_north /= -1) then
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                send_north(j)=work(i,ny-1,k)
            end do
        end do
        call mpi_send(send_north,nx*nz,mpi_real,n_north,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the south
    if(n_south /= -1) then
        call mpi_recv(recv_south,nx*nz,mpi_real,n_south,n_south, &
        pom_comm,istatus,ierr)
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                work(i,1,k)=recv_south(j)
            end do
        end do
    end if

! send ghost cell data to the south
    if(n_south /= -1) then
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                send_south(j)=work(i,2,k)
            end do
        end do
        call mpi_send(send_south,nx*nz,mpi_real,n_south,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the north
    if(n_north /= -1) then
        call mpi_recv(recv_north,nx*nz,mpi_real,n_north,n_north, &
        pom_comm,istatus,ierr)
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                work(i,ny,k)=recv_north(j)
            end do
        end do
    end if

    return
    end subroutine exchange3d_mpi

   


subroutine exchange3d_mpiD(work,nx,ny,nz)
!!  Exchange ghost cells around 3d local grids
!!  one band at a time
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: nx,ny,nz
    real*8 :: work(nx,ny,nz)
    integer :: i,j,k
    integer :: ierr
    integer :: istatus(mpi_status_size)
    real*8 :: send_east(ny*nz),recv_west(ny*nz)
    real*8 :: send_west(ny*nz),recv_east(ny*nz)
    real*8 :: send_north(nx*nz),recv_south(nx*nz)
    real*8 :: send_south(nx*nz),recv_north(nx*nz)

! send ghost cell data to the east
    if(n_east /= -1) then
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                send_east(i)=work(nx-1,j,k)
            end do
        end do
        call mpi_send(send_east,ny*nz,mpi_double_precision,n_east, &
        my_task,pom_comm,ierr)
    end if
! recieve ghost cell data from the west
    if(n_west /= -1) then
        call mpi_recv(recv_west,ny*nz,mpi_double_precision,n_west, &
        n_west,pom_comm,istatus,ierr)
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                work(1,j,k)=recv_west(i)
            end do
        end do
    end if
! send ghost cell data to the west
    if(n_west /= -1) then
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                send_west(i)=work(2,j,k)
            end do
        end do
        call mpi_send(send_west,ny*nz,mpi_double_precision,n_west, &
        my_task,pom_comm,ierr)
    end if
! recieve ghost cell data from the east
    if(n_east /= -1) then
        call mpi_recv(recv_east,ny*nz,mpi_double_precision,n_east, &
        n_east,pom_comm,istatus,ierr)
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                work(nx,j,k)=recv_east(i)
            end do
        end do
    end if
! send ghost cell data to the north
    if(n_north /= -1) then
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                send_north(j)=work(i,ny-1,k)
            end do
        end do
        call mpi_send(send_north,nx*nz,mpi_double_precision,n_north, &
        my_task,pom_comm,ierr)
    end if
! recieve ghost cell data from the south
    if(n_south /= -1) then
        call mpi_recv(recv_south,nx*nz,mpi_double_precision,n_south, &
        n_south,pom_comm,istatus,ierr)
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                work(i,1,k)=recv_south(j)
            end do
        end do
    end if

! send ghost cell data to the south
    if(n_south /= -1) then
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                send_south(j)=work(i,2,k)
            end do
        end do
        call mpi_send(send_south,nx*nz,mpi_double_precision,n_south, &
        my_task,pom_comm,ierr)
    end if
! recieve ghost cell data from the north
    if(n_north /= -1) then
        call mpi_recv(recv_north,nx*nz,mpi_double_precision,n_north, &
        n_north,pom_comm,istatus,ierr)
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                work(i,ny,k)=recv_north(j)
            end do
        end do
    end if
    return
    end subroutine exchange3d_mpiD

    


subroutine order2d_mpi(work2,work4,nx,ny)
!!  Convert a 2nd order 2d matrix to special 4th order 2d matrix
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: nx,ny
    real :: work2(nx,ny),work4(0:nx,0:ny)
    integer :: i,j,k
    integer :: ierr
    integer :: istatus(mpi_status_size)
    real :: send_east(ny),recv_west(ny)
    real :: send_north(nx),recv_south(nx)

    work4=0.
    do i=1,nx
        do j=1,ny
            work4(i,j)=work2(i,j)
        end do
    end do

! send ghost cell data to the east
    if(n_east /= -1) then
        do j=1,ny
            send_east(j)=work2(nx-2,j)
        end do
        call mpi_send(send_east,ny,mpi_real,n_east,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the west
    if(n_west /= -1) then
        call mpi_recv(recv_west,ny,mpi_real,n_west,n_west, &
        pom_comm,istatus,ierr)
        do j=1,ny
            work4(0,j)=recv_west(j)
        end do
    end if

! send ghost cell data to the north
    if(n_north /= -1) then
        do i=1,nx
            send_north(i)=work2(i,ny-2)
        end do
        call mpi_send(send_north,nx,mpi_real,n_north,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the south
    if(n_south /= -1) then
        call mpi_recv(recv_south,nx,mpi_real,n_south,n_south, &
        pom_comm,istatus,ierr)
        do i=1,nx
            work4(i,0)=recv_south(i)
        end do
    end if

    return
    end subroutine order2d_mpi

    


subroutine order3d_mpi(work2,work4,nx,ny,nz)
!!  Convert a 2nd order 3d matrix to special 4th order 3d matrix
    implicit none
    include 'mpif.h'
    include 'pom.h'
    integer :: nx,ny,nz
    real :: work2(nx,ny,nz),work4(0:nx,0:ny,nz)
    integer :: i,j,k
    integer :: ierr
    integer :: istatus(mpi_status_size)
    real :: send_east(ny*nz),recv_west(ny*nz)
    real :: send_north(nx*nz),recv_south(nx*nz)

    work4=0.
    do i=1,nx
        do j=1,ny
            do k=1,nz
                work4(i,j,k)=work2(i,j,k)
            end do
        end do
    end do

! send ghost cell data to the east
    if(n_east /= -1) then
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                send_east(i)=work2(nx-2,j,k)
            end do
        end do
        call mpi_send(send_east,ny*nz,mpi_real,n_east,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the west
    if(n_west /= -1) then
        call mpi_recv(recv_west,ny*nz,mpi_real,n_west,n_west, &
        pom_comm,istatus,ierr)
        do k=1,nz
            do j=1,ny
                i=j+(k-1)*ny
                work4(0,j,k)=recv_west(i)
            end do
        end do
    end if

! send ghost cell data to the north
    if(n_north /= -1) then
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                send_north(j)=work2(i,ny-2,k)
            end do
        end do
        call mpi_send(send_north,nx*nz,mpi_real,n_north,my_task, &
        pom_comm,ierr)
    end if
! recieve ghost cell data from the south
    if(n_south /= -1) then
        call mpi_recv(recv_south,nx*nz,mpi_real,n_south,n_south, &
        pom_comm,istatus,ierr)
        do k=1,nz
            do i=1,nx
                j=i+(k-1)*nx
                work4(i,0,k)=recv_south(j)
            end do
        end do
    end if

    return
    end subroutine order3d_mpi

