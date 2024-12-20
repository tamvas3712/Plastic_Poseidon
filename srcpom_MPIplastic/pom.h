! pom.h

! contain parameters for the model domain and for the decomposition
! local domain size, and common POM variables

    include 'ibmall.fcm'
! ______________________________________________________________________
! Grid size parameters
    integer :: &
    im_global      ,& !  number of global grid points in x
    jm_global      ,& !  number of global grid points in y
    kb             ,& !  number of grid points in z
    im_local       ,& !  number of local grid points in x
    jm_local       ,& !  number of local grid points in y
    n_proc         ,& !  number of processors
    nsi_global     ,& !  number of global nsi !changing
    nsi_local      ,& ! number of local nsi  !changing
    n_west         ,& !  western parallel processor ID
    n_east         ,& !  eastern parallel processor ID
    n_south        ,& !  southern parallel processor ID
    n_north         ! northern parallel processor ID


! Correct values for im_local and jm_local are found using
!   n_proc=(im_global-2)/(im_local-2)*(jm_global-2)/(jm_local-2)
! Values higher than necessary will not cause the code to fail, but
! will allocate more memory than is necessary. Value that are too low
! will cause the code to exit


    parameter( &
    im_global=431  , &
    jm_global=156   , &
    kb=25          , &

    im_local=431    , &
    jm_local=156    , &

    nsi_global =500000, &
    nsi_local = 100000, &
    n_proc=5       )

! ______________________________________________________________________
! Efective grid size
    integer :: &
    im             ,& !  number of grid points used in local x domains
    imm1           ,& !  im-1
    imm2           ,& !  im-2
    jm             ,& !  number of grid points used in local y domains
    jmm1           ,& !  jm-1
    jmm2           ,& !  jm-2
    kbm1           ,& !  kb-1
    kbm2            ! kb-2

! Note that im and jm may be different between local domains
! im and jm are equal to or lower than im_local and jm_local, depending
! on the use of correct values for im_local and jm_local

    common/blksiz/ &
    im             , &
    imm1           , &
    imm2           , &
    jm             , &
    jmm1           , &
    jmm2           , &
    kbm1           , &
    kbm2

! ______________________________________________________________________
! Parallel variables
    integer :: &
    my_task        ,& !  actual parallel processor ID
    master_task    ,& !  master processor ID
    pom_comm       ,& !  POM model MPI group communicator
    i_global       ,& !  global i index for each point in local domain
    j_global       ,& !  global j index for each point in local domain
    n_global        ! global n index for IBM KTSIARAS


    common/blkpar/ &
    my_task        , &
    master_task    , &
    pom_comm       , &
    i_global(im_local), &
    j_global(jm_local), &
    n_global(nmax) , &
    n_west         , &
    n_east         , &
    n_south        , &
    n_north



! ______________________________________________________________________
! Scalars
    integer :: &
    iint           , &
    iprint         ,& !  interval in iint at which variables are printed
    iend           ,& !  total internal mode time steps
    iswtch         ,& !  time step interval to switch from prtd1 to prtd2
    irestart       , &
    error_status   , &
    iyr,iyr1,imon,iday   !KTSIARAS

    real :: &
    dti            ,& !  internal (3-D) time step (s)
    dti2           ,& !  2*dti
    grav           ,& !  gravity constant (S.I. units)
    pi             ,& !  pi
    rhoref         ,& !  reference density
    sbias          ,& !  salinity bias
    small          ,& !  small value
    tbias          ,& !  temperature bias
    time           ,& !  model time (days)
    write_rst      , &
    days           ,& !  run duration in days
    horcon         ,& !  smagorinsky diffusivity coefficient
    prtd1          ,& !  initial print interval (days)
    prtd2          ,& !  final print interval (days)
    swtch          ,& !  time to switch from prtd1 to prtd2
    time0             !  initial time (days)

    common/blkcon/ &
    dti            , &
    dti2           , &
    grav           , &
    pi             , &
    rhoref         , &
    sbias          , &
    small          , &
    tbias          , &
    time           , &
    write_rst      , &
    iint           , &
    iprint         , &
    days           , &
    horcon         , &
    prtd1          , &
    prtd2          , &
    swtch          , &
    time0          , &
    iend           , &
    iswtch         , &
    irestart       , &
    error_status   , &
    iyr,imon,iday   !KTSIARAS


! ______________________________________________________________________
! 1-D arrays
    real :: &
    dz             ,& !  z(k)-z(k+1)
    dzz            ,& !  zz(k)-zz(k+1)
    z              ,& !  sigma coordinate from z=0 (surface) to z=-1 (bottom)
    zz              ! sigma coordinate, intermediate between z

    common/blk1d/ &
    dz(kb)         , &
    dzz(kb)        , &
    z(kb)          , &
    zz(kb)

! ______________________________________________________________________
! 2-D arrays
    real :: &
    art            ,& !  cell area centered on T grid points
    aru            ,& !  cell area centered on U grid points
    arv            ,& !  cell area centered on V grid points
    dum            ,& !  mask for u velocity
    dvm            ,& !  mask for v velocity
    dx             ,& !  grid spacing in x
    dy             ,& !  grid spacing in y
    east_e         ,& !  horizontal coordinate of elevation points in x
    fsm            ,& !  mask for scalar variables
    h              ,& !  bottom depth
    dt             ,& !  =h+ssh
    north_e        ,& !  horizontal coordinate of elevation points in y
    wubot          ,& !  x-momentum flux at the bottom
    wvbot          ,& !  y-momentum flux at the bottom
    cbc            ,& !  bottom friction coefficient
!--KTSIARAS atmos for bulk folmulae
    uairf          , &
    vairf          , &
    stokesxf       , &
    stokesyf       , &
    waveperf       , &
    waveheif       , &
    uairb          , &
    vairb          , &
    uair           , &
    vair           , &
    stokesx        , &
    stokesy        , &
    waveper        , &
    wavehei        , &
    alon           , &
    alat           , &
    rbc
    integer :: lbc               !boundary indices for particles reflections
    integer :: iw2d


    common/blk2d/ &
    cbc(im_local,jm_local)     , &
    dum(im_local,jm_local)     , &
    dvm(im_local,jm_local)     , &
    dx(im_local,jm_local)      , &
    dy(im_local,jm_local)      , &
    east_e(im_local,jm_local)  , &
    fsm(im_local,jm_local)     , &
    h(im_local,jm_local)       , &
    dt(im_local,jm_local)       , &
    north_e(im_local,jm_local) , &
    wubot(im_local,jm_local)   , &
    wvbot(im_local,jm_local)   , &
!--KTSIARAS atmos for bulk folmulae
    uairf(im_local,jm_local,9)   , &
    vairf(im_local,jm_local,9)   , &
    stokesxf(im_local,jm_local,9), &
    stokesyf(im_local,jm_local,9), &
    waveperf(im_local,jm_local,9), &
    waveheif(im_local,jm_local,9), &
    uair(im_local,jm_local)    , &
    vair(im_local,jm_local)    , &
    stokesx(im_local,jm_local) , &
    stokesy(im_local,jm_local) , &
    waveper(im_local,jm_local) , &
    wavehei(im_local,jm_local) , &
    alon(im_local,jm_local)     , &
    alat(im_local,jm_local)     , &
    iw2d(im_local,jm_local)     , &
    lbc(im_local,jm_local), &
    rbc(im_local,jm_local)

! ______________________________________________________________________
! 3-D arrays
    real :: &
    aam            ,& !  horizontal kinematic viscosity
    kh             ,& !  vertical diffusivity
    rho            ,& !  density
    rmean          ,& !  horizontally averaged density
    elb            ,& !  sea surface elevation
    s              ,& !  salinity at time n
    t              ,& !  temperature at time n
    u              ,& !  horizontal velocity in x at time n
    v              ,& !  horizontal velocity in y at time n
    w              ,& !  sigma coordinate vertical velocity
!--KTSIARAS-------------------------------------------------------------
    bac3d             !  bacteria

    common/blk3d/ &
    aam(im_local,jm_local,kb)  , &
    kh(im_local,jm_local,kb)   , &
    rho(im_local,jm_local,kb)  , &
    rmean(im_local,jm_local,kb), &
    elb(im_local,jm_local), &
    s(im_local,jm_local,kb)    , &
    t(im_local,jm_local,kb)    , &
    u(im_local,jm_local,kb)    , &
    v(im_local,jm_local,kb)    , &
    w(im_local,jm_local,kb)    , &
!--KTSIARAS-------------------------------------------------------------
    bac3d(im_local,jm_local,kb-1)

! ______________________________________________________________________
! Character variables
    character(26) :: &
    time_start      ! date and time of start of initial run of model

    character(40) :: &
    source         , &
    title

    character(120) :: &
    netcdf_file    , &
    read_rst_file  , &
    write_rst_file

    common/blkchar/ &
    time_start     , &
    source         , &
    title          , &
    netcdf_file    , &
    read_rst_file  , &
    write_rst_file

