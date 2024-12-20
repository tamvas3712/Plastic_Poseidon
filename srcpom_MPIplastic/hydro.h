    integer :: nfre,nfre_eco,nt,nt_eco
    parameter (nfre=3)                 !daily frequency of Hydrodynamic output
    PARAMETER (NT=24/nfre+1)  !number of daily BCs
!------------------------------------------------------------

    real :: &
    u_f          ,& !  u-velolity
    v_f          ,& !  v-velolity
    w_f          ,& !  w-velolity
    t_f          ,& !  temperature
    s_f          ,& !  salinity
    aam_f        ,& !  horizontal diffusion
    kh_f         ,& !  vertical diffusivity
    el_f          ! sea surface elevation


    common/hydro/ &
    u_f(im_local,jm_local,kb,nt), &
    v_f(im_local,jm_local,kb,nt), &
    w_f(im_local,jm_local,kb,nt), &
    t_f(im_local,jm_local,kb,nt), &
    s_f(im_local,jm_local,kb,nt), &
    aam_f(im_local,jm_local,kb,nt), &
    kh_f(im_local,jm_local,kb,nt), &
    el_f(im_local,jm_local,nt)
