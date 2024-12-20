! io_pnetcdf.F

! parallel NetCDF I/O


subroutine def_var_pnetcdf(ncid,name,nvdims,vdims,varid, &   !netcdf variable definition (general)
    long_name,units,coords,lcoords)
!! -?? need to define
    # include "mpif.h"
    # include "pnetcdf.inc"
    integer :: vdims(nvdims)
    integer :: ncid,nvdims,varid
    logical :: lcoords
    character*(*) name,long_name,units,coords
    integer :: status
    integer(mpi_offset_kind) length
! define variable
    status=nfmpi_def_var(ncid,name,nf_float,nvdims,vdims,varid)
    call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)
! define attributes
    length=len(trim(long_name))
    status=nfmpi_put_att_text(ncid,varid,'long_name', &
    length,trim(long_name))
    call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)

    length=len(trim(units))
    status=nfmpi_put_att_text(ncid,varid,'units',length,trim(units))
    call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)
! add coordinates attribute, if necessary
    if(lcoords) then
        length=len(trim(coords))
        status=nfmpi_put_att_text(ncid,varid,'coordinates', &
        length,trim(coords))
        call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)
    end if
    return
    end subroutine def_var_pnetcdf

    
subroutine handle_error_pnetcdf(routine,status,nf_noerr) !netcdf error (general)
!! Subroutine to handle error.
    implicit none
    include 'pom.h'
    integer :: status,nf_noerr
    character*(*) routine
    if(status /= nf_noerr) then
        error_status=1
        if(my_task == master_task) write(*,'(/a,a)') &
        'Error: NetCDF routine ',routine
    end if
    return
    end subroutine handle_error_pnetcdf


subroutine write_ibmrestart_pnetcdf !KTSIARAS writes ibm daily restart file (out/pom.ibmplastic.nc)
!!  Writes ibm 1-d array
!!  Write data for a seamless restart
    implicit none
    # include "mpif.h"
    # include "pnetcdf.inc"
    include 'pom.h'
    character(120) :: netcdf_out_file,str_tmp
    integer :: nprint
    integer :: ncid,status
    integer :: nsi_dimid,species_id,xfi_id,nfi_id,variant_id,xsi_id,ysi_id, &
    iflagVariant_id,weirate_id,rbufcum_id, &
    conxfi_id,ref_id,iarea_id,indsi_id, &
    icountEgg_id,keepeggs_id,adday_coho_id, &
    icountdays_id,fdep_id,rbt_id,fou_id,nschoo_id,dista_id

    integer :: vdims(1)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start(1),edge(1)
    integer :: n,k,i
    real :: totmacro1(6,2)

    write(*,*)'Start writing ibm'
! create netcdf restart file
!      nprint=(iint+time0*86400/dti)/irestart
    nprint=iint/irestart
    write(netcdf_out_file,'(''../out/'',a,''.'',''ibmplastic.nc'')') &
    trim(write_rst_file)
    if(my_task == 0) write(*,'(/''writing file '',a)') &
    trim(netcdf_out_file)
    status=nfmpi_create(pom_comm,netcdf_out_file, &
    nf_clobber+nf_64bit_offset,mpi_info_null,ncid)
    call handle_error_pnetcdf('nf_create',status,nf_noerr)

! define global attributes
    length=len(trim(title))
    status=nfmpi_put_att_text(ncid,nf_global,'title',length, &
    trim(title))
    call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)

    str_tmp='ibm restart file'
    length=len(trim(str_tmp))
    status=nfmpi_put_att_text(ncid,nf_global,'description',length, &
    trim(str_tmp))
    call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)

! collect nsi from all CPUs new_nsi_global=sum(nsi)
!      call MPI_GATHER(SENDBUF, n_proc, MPI_INT, RECVBUF, 1,
!        MPI_INT, 0, COMM, IERROR)

! define dimensions
    length=nsi_global
    status=nfmpi_def_dim(ncid,'nsi',length,nsi_dimid)
    call handle_error_pnetcdf('nf_def_dim',status,nf_noerr)

! define variables and their attributes
    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'species',1,vdims,species_id, &
    'species', &
    ' ', &
    'SIs',.true.)

!      vdims(1)=nsi_dimid
!      call def_var_pnetcdf(ncid,'xfi',1,vdims,xfi_id,
!     $                     'weight',
!     $                     'gr',
!     $                     'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'nfi',1,vdims,nfi_id, &
    'number of individuals', &
    ' ', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'variant',1,vdims,variant_id, &
    'variant', &
    ' ', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'xsi',1,vdims,xsi_id, &
    'Longitude', &
    'deg. East', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'ysi',1,vdims,ysi_id, &
    'Latitude', &
    'deg. North', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'ref',1,vdims,ref_id, &
    'flag indicating status of SI', &
    'ref=-2 outside boundary, ref>0 land etc', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'iarea',1,vdims,iarea_id, &
    'flag indicating area of origin', &
    'flag', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'icountdays',1,vdims,icountdays_id, &
    'count of days at sea', &
    '#days', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'fdep',1,vdims,fdep_id, &
    'depth in the water column', &
    'm', &
    'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'rbt',1,vdims,rbt_id, &
    'hours spent on beach', &
    'h', &
    'SI retention time',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'fou',1,vdims,fou_id, &
    'biofouling thickness', &
    'm', &
    'SI location',.true.)

!      vdims(1)=nsi_dimid
!      call def_var_pnetcdf(ncid,'nschoo',1,vdims,nschoo_id,
!     $                     '#of gathered SIs',
!     $                     '#',
!     $                     'SI location',.true.)

    vdims(1)=nsi_dimid
    call def_var_pnetcdf(ncid,'dista',1,vdims,dista_id, &
    'travelled distance', &
    'm', &
    'SI location',.true.)

! end definitions
    status=nfmpi_enddef(ncid)
    call handle_error_pnetcdf('nf_enddef',status,nf_noerr)
      

! write data
    start(1)=n_global(1)
    edge(1)=nsi
    status=nfmpi_put_vara_int_all(ncid,species_id,start,edge,species)
    call handle_error_pnetcdf('nfmpi_put_vara_int_all', &
    status,nf_noerr)
!      status=nfmpi_put_vara_real_all(ncid,xfi_id,start,edge,xfi)
!      call handle_error_pnetcdf('nfmpi_put_vara_real_all',
!     $                                                  status,nf_noerr)
    status=nfmpi_put_vara_real_all(ncid,nfi_id,start,edge,nfi)
    call handle_error_pnetcdf('nfmpi_put_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_put_vara_int_all(ncid,variant_id,start,edge,variant)
    call handle_error_pnetcdf('nfmpi_put_vara_int_all', &
    status,nf_noerr)
    status=nfmpi_put_vara_real_all(ncid,xsi_id,start,edge,xsi)
    call handle_error_pnetcdf('nfmpi_put_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_put_vara_real_all(ncid,ysi_id,start,edge,ysi)
    call handle_error_pnetcdf('nfmpi_put_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_put_vara_int_all(ncid,ref_id,start,edge,ref)
    call handle_error_pnetcdf('nfmpi_put_vara_int_all', &
    status,nf_noerr)
    status=nfmpi_put_vara_int_all(ncid,iarea_id,start,edge,iarea)
    call handle_error_pnetcdf('nfmpi_put_vara_int_all', &
    status,nf_noerr)

!      status=nfmpi_put_vara_int_all(ncid,nschoo_id,start,edge,nschoo)
!      call handle_error_pnetcdf('nfmpi_put_vara_int_all',
!     $                                                  status,nf_noerr)

!--update days at sea
    do k = 1, nsi
        if(variant(k) > 0)then
            icountdays(k)=icountdays(k)+1
        endif
    enddo

    status=nfmpi_put_vara_int_all(ncid,icountdays_id,start,edge, &
    icountdays)
    call handle_error_pnetcdf('nfmpi_put_vara_int_all', &
    status,nf_noerr)
    status=nfmpi_put_vara_real_all(ncid,fdep_id,start,edge,fdep)
    call handle_error_pnetcdf('nfmpi_put_vara_real_all', &
    status,nf_noerr)

    status=nfmpi_put_vara_real_all(ncid,rbt_id,start,edge,rbt)
    call handle_error_pnetcdf('nfmpi_put_vara_real_all', &
    status,nf_noerr)

    status=nfmpi_put_vara_real_all(ncid,fou_id,start,edge,fou)
    call handle_error_pnetcdf('nfmpi_put_vara_real_all', &
    status,nf_noerr)

    status=nfmpi_put_vara_real_all(ncid,dista_id,start,edge,dista)
    call handle_error_pnetcdf('nfmpi_put_vara_real_all', &
    status,nf_noerr)



! close file
    status=nfmpi_close(ncid)
    call handle_error_pnetcdf('nf_close',status,nf_noerr)

    do n=1,2
        do i=1,6
            totmacro1(i,n)=0.
        enddo
    enddo

    do k=1,nsi
        i=variant(k)
        n=species(k)
        if(variant(k) /= 0) &
        totmacro1(i,n)=totmacro1(i,n)+nfi(k)
    enddo

    do n=1,2
        do i=1,6
            write(*,*)'TOTMACRO1',my_task,n,i,totmacro1(i,n)
        enddo
    enddo


    return
    end subroutine write_ibmrestart_pnetcdf


subroutine read_grid_pnetcdf  !reads model grid, bathymetry etc (pom.grid.nc)
    implicit none
    # include "mpif.h"
    # include "pnetcdf.inc"
    include 'pom.h'
    character(120) :: netcdf_grid_file
    integer :: z_varid,zz_varid,dx_varid,dy_varid,east_c_varid, &
    east_e_varid,east_u_varid,east_v_varid,north_c_varid, &
    north_e_varid,north_u_varid,north_v_varid, &
    h_varid,fsm_varid,dum_varid,dvm_varid
    integer :: ncid,status,i,j
    integer :: vdims(2)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start(2),edge(2)

! open netcdf file
    write(netcdf_grid_file,'(''../in/'',a,''.grid.nc'')') &
    trim(netcdf_file)

    if(my_task == 0) write(*,'(/''reading file POMplastic Grid'',a)') &
    trim(netcdf_grid_file)
    status=nfmpi_open(pom_comm,netcdf_grid_file,nf_nowrite, &
    mpi_info_null,ncid)
    call handle_error_pnetcdf(netcdf_grid_file, &
    status,nf_noerr)

! get variables
    status=nfmpi_inq_varid(ncid,'z',z_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: z',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'zz',zz_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: zz',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'dx',dx_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: dx',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'dy',dy_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: dy',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'east_e',east_e_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: east_e',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'north_e',north_e_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: north_e',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'h',h_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: h',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'fsm',fsm_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: fsm',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'dum',dum_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: dum',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'dvm',dvm_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: dvm',status,nf_noerr)

! get data
    start(1)=1
    edge(1)=kb
    status=nfmpi_get_vara_real_all(ncid,z_varid,start,edge,z)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,zz_varid,start,edge,zz)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)

    start(1)=i_global(1)
    start(2)=j_global(1)
    edge(1)=im
    edge(2)=jm
    status=nfmpi_get_vara_real_all(ncid,dx_varid,start,edge,dx)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,dy_varid,start,edge,dy)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,east_e_varid,start,edge,east_e)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,north_e_varid,start,edge,north_e)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,h_varid,start,edge,h)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,fsm_varid,start,edge,fsm)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,dum_varid,start,edge,dum)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,dvm_varid,start,edge,dvm)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all',status,nf_noerr)
! close file:
    status=nfmpi_close(ncid)
    call handle_error_pnetcdf('nf_close: grid',status,nf_noerr)

    return
    end subroutine read_grid_pnetcdf




subroutine read_ibm_pnetcdf !KTSIARAS reads ibm restart file (in/pom.ibmplastic.nc)
!!  Read ibm attributes
    implicit none
    # include "mpif.h"
    # include "pnetcdf.inc"
    include 'pom.h'

    character(120) :: netcdf_ic_file
    integer :: ncid,status
    integer :: nsi_id,species_id,xfi_id,nfi_id,cfi_id,variant_id,xsi_id,ysi_id, &
    iflagVariant_id,weirate_id,zoop1_id,zoop2_id,zoop3_id,rbufcum_id, &
    temp_id,grdot_id,conxfi_id,ref_id,iarea_id,indsi_id,topo_id, &
    icountEgg_id,keepeggs_id,adday_coho_id,speed_id,fdep_id, &
    icountdays_id,rbt_id,fou_id,nschoo_id

    integer :: vdims(4)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start(1),edge(1)
    integer :: is,k,i,j,n
    real :: totmacro1(6,2)

! open netcdf file
    write(netcdf_ic_file,'(''../in/'',a,''.ibmplastic.nc'')') &
    trim(netcdf_file)
    if(my_task == 0) write(*,'(/''reading file '',a)') &
    trim(netcdf_ic_file)
    status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite, &
    mpi_info_null,ncid)
    call handle_error_pnetcdf(netcdf_ic_file, &
    status,nf_noerr)

! get variables
    status=nfmpi_inq_varid(ncid,'species',species_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:spe',status,nf_noerr)
!      status=nfmpi_inq_varid(ncid,'xfi',xfi_id)
!      call handle_error_pnetcdf('nfmpi_inq_varid:xfi',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'nfi',nfi_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:nfi',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'variant',variant_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:variant',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'xsi',xsi_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:xsi',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'ysi',ysi_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:ysi',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'ref',ref_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:ref',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'iarea',iarea_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:iarea',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'fdep',fdep_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:fdep',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'icountdays',icountdays_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:icountdays',status, &
    nf_noerr)
    status=nfmpi_inq_varid(ncid,'rbt',rbt_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:rbt',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'fou',fou_id)
    call handle_error_pnetcdf('nfmpi_inq_varid:fou',status,nf_noerr)
!      status=nfmpi_inq_varid(ncid,'nschoo',nschoo_id)
!      call handle_error_pnetcdf('nfmpi_inq_varid:nschoo',
!     &status,nf_noerr)



! get data
    start(1)=n_global(1)
    edge(1)=nsi
    status=nfmpi_get_vara_int_all(ncid,species_id,start,edge,species)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
!      status=nfmpi_get_vara_real_all(ncid,xfi_id,start,edge,xfi)
!      call handle_error_pnetcdf('nfmpi_get_vara_real_all',
!     $                                                  status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,nfi_id,start,edge,nfi)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_int_all(ncid,variant_id,start,edge,variant)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,xsi_id,start,edge,xsi)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,ysi_id,start,edge,ysi)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_int_all(ncid,ref_id,start,edge,ref)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_int_all(ncid,iarea_id,start,edge,iarea)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_int_all(ncid,icountdays_id,start,edge, &
    icountdays)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,rbt_id,start,edge,rbt)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,fdep_id,start,edge,fdep)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,fou_id,start,edge,fou)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
!      status=nfmpi_get_vara_int_all(ncid,nschoo_id,start,edge,nschoo)
!      call handle_error_pnetcdf('nfmpi_get_vara_real_all',
!     $                                                  status,nf_noerr)


! close file
    status=nfmpi_close(ncid)
    call handle_error_pnetcdf('nf_close',status,nf_noerr)

    do n=1,2
        do k=1,6
            totmacro1(k,n)=0.
        enddo
    enddo

    do k=1,nsi
        i=variant(k)
        n=species(k)
        if(variant(k) /= 0) &
        totmacro1(i,n)=totmacro1(i,n)+nfi(k)
    enddo

    do n=1,2
        do i=1,6
            write(*,*)'TOTMACRO0',my_task,n,i,totmacro1(i,n)
        enddo
    enddo


!      do k=1,nsi
!      if(species(k).ne.0..and.iarea(k).gt.1)then
!      write(*,*)'READIBM',species(k),variant(k),nfi(k),ref(k),iarea(k)
!      endif
!      enddo
!      write(*,*)'TEST-IBMREAD',my_task,xfi(1),variant(1),xsi(1)

    return
    end subroutine read_ibm_pnetcdf


subroutine read_bac3d_pnetcdf  !KTSIARAS reads bacteria 3-d field (in/pom.ecologyibm.nc), used for biofouling
!!  Read ecology 3-d state (includes land)
    implicit none
    # include "mpif.h"
    # include "pnetcdf.inc"
    include 'pom.h'
    character(120) :: netcdf_ic_file
    integer :: ncid,status
    integer :: varid,varbid
    integer :: vdims(4)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start1(4),edge1(4),start2(3),edge2(3)
    integer :: is,k,i,j
! open netcdf file
    write(netcdf_ic_file,'(''../in/'',a,''.ecologyibm.nc'')') &
    trim(netcdf_file)
    if(my_task == 0) write(*,'(/''reading file '',a)') &
    trim(netcdf_ic_file)
    status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite, &
    mpi_info_null,ncid)
    call handle_error_pnetcdf(netcdf_ic_file, &
    status,nf_noerr)

! get variables
    status=nfmpi_inq_varid(ncid,'bact',varid)
    call handle_error_pnetcdf('nfmpi_inq_varid:bact',status,nf_noerr)

! get data

    start2(1)=i_global(1)
    start2(2)=j_global(1)
    start2(3)=1
    edge2(1)=im
    edge2(2)=jm
    edge2(3)=kb-1
    status=nfmpi_get_vara_real_all(ncid,varid,start2,edge2,bac3d)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all-bact', &
    status,nf_noerr)

!      write(*,*)im,jm,kb-1,p_state,i_global(1),j_global(1)
!      do is=1,p_state
!      do k=1,kb-1
!       write(*,*)'TESTREAD',is,k,varp(340,40,k,is)
!      enddo
!      enddo

! close file
    status=nfmpi_close(ncid)
    call handle_error_pnetcdf('nf_close',status,nf_noerr)

    return
    end subroutine read_bac3d_pnetcdf


subroutine read_hydroibm_pnetcdf!KTSIARAS reads 3-hour hydro fields (elevation, U,V,W,T,S,KH) (in/pom.hydroibm.nc)
    implicit none
    # include "mpif.h"
    # include "pnetcdf.inc"
    include 'pom.h'
    include 'hydro.h'
    character(120) :: netcdf_ic_file
    integer :: ncid,status,n
    integer :: varid,varbid
    integer :: vdims(3)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start(4),edge(4),start2(3),edge2(3)

! open netcdf file
    write(netcdf_ic_file,'(''../in/'',a,''.hydroibm.nc'')') &
    trim(netcdf_file)
    if(my_task == 0) write(*,'(/''reading file '',a)') &
    trim(netcdf_ic_file)
    status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite, &
    mpi_info_null,ncid)
    call handle_error_pnetcdf(netcdf_ic_file, &
    status,nf_noerr)

! get variables

    write(*,*)'hydroibm',my_task

    status=nfmpi_inq_varid(ncid,'u',varid)
    call handle_error_pnetcdf('nfmpi_inq_varid:u',status,nf_noerr)

    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    start(4)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=kb
    edge(4)=nt
    status=nfmpi_get_vara_real_all(ncid,varid,start,edge,u_f)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)


    status=nfmpi_inq_varid(ncid,'v',varid)
    call handle_error_pnetcdf('nfmpi_inq_varid:v_f',status,nf_noerr)

    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    start(4)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=kb
    edge(4)=nt
    status=nfmpi_get_vara_real_all(ncid,varid,start,edge,v_f)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)


    status=nfmpi_inq_varid(ncid,'w',varid)
    call handle_error_pnetcdf('nfmpi_inq_varid:w_f',status,nf_noerr)

    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    start(4)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=kb
    edge(4)=nt
    status=nfmpi_get_vara_real_all(ncid,varid,start,edge,w_f)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)


    status=nfmpi_inq_varid(ncid,'t',varid)
    call handle_error_pnetcdf('nfmpi_inq_varid:t_f',status,nf_noerr)

    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    start(4)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=kb
    edge(4)=nt
    status=nfmpi_get_vara_real_all(ncid,varid,start,edge,t_f)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)

    write(*,*)'hydroibm-t',my_task

    status=nfmpi_inq_varid(ncid,'kh',varid)
    call handle_error_pnetcdf('nfmpi_inq_varid:kh_f',status,nf_noerr)

    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    start(4)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=kb
    edge(4)=nt
    status=nfmpi_get_vara_real_all(ncid,varid,start,edge,kh_f)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)

    write(*,*)'hydroibm-kh',my_task


    status=nfmpi_inq_varid(ncid,'aam',varid)
    call handle_error_pnetcdf('nfmpi_inq_varid:aam_f',status,nf_noerr)

    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    start(4)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=kb
    edge(4)=nt
    status=nfmpi_get_vara_real_all(ncid,varid,start,edge,aam_f)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)

    write(*,*)'hydroibm-aam',my_task

! close file
    status=nfmpi_close(ncid)
    call handle_error_pnetcdf('nf_close',status,nf_noerr)

    do n=1,nt
        write(*,*)'READ-HYDRO',my_task,n,t_f(340,150,1,n), &
        u_f(340,150,1,n),v_f(340,150,1,n),kh_f(340,150,2,n), &
        aam_f(340,150,1,n)
    enddo

    return
    end subroutine read_hydroibm_pnetcdf


subroutine read_wave_pnetcdf(n) !KTSIARAS reads 3-hour wave fields (stokes drift, wave height/period)
!!  read wind stress
    implicit none
    # include "mpif.h"
    # include "pnetcdf.inc"
    include 'pom.h'
    integer :: n
    integer :: i,j
    character(120) :: netcdf_file2
    integer :: ncid,status
    integer :: waveper_varid,wavehei_varid,wavedir_varid, &
    stokesx_varid,stokesy_varid
    integer :: vdims(2)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start(3),edge(3)

! open netcdf file
    write(netcdf_file2,'(''../in/'',a,''.wave.nc'')')trim(netcdf_file)

    if(my_task == 0) write(*,'(/''reading file '',a)') &
    trim(netcdf_file2)
    write(*,*)'TEST-n',n
    status=nfmpi_open(pom_comm,netcdf_file2,nf_nowrite, &
    mpi_info_null,ncid)
    call handle_error_pnetcdf(netcdf_file2, &
    status,nf_noerr)

! get variables
    status=nfmpi_inq_varid(ncid,'waveper',waveper_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: waveper', &
    status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'wavehei',wavehei_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: wavehei', &
    status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'stokesx',stokesx_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: stokesx', &
    status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'stokesy',stokesy_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: stokesy', &
    status,nf_noerr)




! get data
    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=9
    status=nfmpi_get_vara_real_all(ncid,waveper_varid,start,edge, &
    waveperf)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,wavehei_varid,start,edge, &
    waveheif)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,stokesx_varid,start,edge, &
    stokesxf)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,stokesy_varid,start,edge, &
    stokesyf)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)


! close file
    status=nfmpi_close(ncid)
    call handle_error_pnetcdf('nf_close',status,nf_noerr)

    return
    end subroutine read_wave_pnetcdf

    subroutine read_wind2_pnetcdf(n) !KTSIARAS reads 3-hour wind fields
! read wind stress
    implicit none
    # include "mpif.h"
    # include "pnetcdf.inc"
    include 'pom.h'
    integer :: n
    integer :: i,j
    character(120) :: netcdf_file2
    integer :: ncid,status
    integer :: uair_varid,vair_varid,tair_varid,rhum_varid,sol_varid, &
    rlond_varid,rain_varid
    integer :: vdims(2)
    integer(mpi_offset_kind) length
    integer(mpi_offset_kind) start(3),edge(3)

! open netcdf file
!      write(netcdf_file2,'(''../in/'',a,''.atmos.'',i4.4,''.nc'')')
!     $ trim(netcdf_file),n
    write(netcdf_file2,'(''../in/'',a,''.atmos.nc'')')trim(netcdf_file)

    if(my_task == 0) write(*,'(/''reading file '',a)') &
    trim(netcdf_file2)
    write(*,*)'TEST-n',n
    status=nfmpi_open(pom_comm,netcdf_file2,nf_nowrite, &
    mpi_info_null,ncid)
    call handle_error_pnetcdf(netcdf_file2, &
    status,nf_noerr)

! get variables
    status=nfmpi_inq_varid(ncid,'uair',uair_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: uair',status,nf_noerr)
    status=nfmpi_inq_varid(ncid,'vair',vair_varid)
    call handle_error_pnetcdf('nfmpi_inq_varid: vair',status,nf_noerr)

! get data
    start(1)=i_global(1)
    start(2)=j_global(1)
    start(3)=1
    edge(1)=im
    edge(2)=jm
    edge(3)=9
    status=nfmpi_get_vara_real_all(ncid,uair_varid,start,edge,uairf)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)
    status=nfmpi_get_vara_real_all(ncid,vair_varid,start,edge,vairf)
    call handle_error_pnetcdf('nfmpi_get_vara_real_all', &
    status,nf_noerr)

! close file
    status=nfmpi_close(ncid)
    call handle_error_pnetcdf('nf_close',status,nf_noerr)

    return
    end subroutine read_wind2_pnetcdf

