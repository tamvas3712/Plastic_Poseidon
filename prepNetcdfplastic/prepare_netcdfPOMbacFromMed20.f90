program prepare_netcdf_bac
!!  This program reads INIT.DAT and bac3d and writes netcdf data ??
    include "netcdf.inc"
    include 'grid.h'
    include 'sesame.fcm'
    include 'mappings.h'
    include 'fsm3d.h'
    character date*6
    real :: VAR(IM,JM,KB-1,P_STATE),VARB(IM,JM,B_STATE)
    real :: var1d(N_COMPg),var1db(n_upper$g)
    real :: VAR3D(IM,JM,KB-1),VAR2D(IM,JM)
    real :: var1(im-1,jm-1),elb1(im-1,jm-1)
    real :: bac3d(IM,JM,KB-1),chl3d(IM,JM,KB-1)
    real :: bac3dav(IM,JM,KB-1),chl3dav(IM,JM,KB-1)
    real :: bac3dmo(IM,JM,KB-1,12),chl3dmo(IM,JM,KB-1,12)
    real :: bac3dz(im,jm,klev),chl3dz(im,jm,klev)
    parameter (MINCHK=-900)
    integer :: NCID1,KLEVID,JMID,IMID,PSTATEID,BSTATEID
    character CFILE*128,DATENAME*6,indf1*2
    integer :: ldf1

!---read RESTART
     
!---read INIT.DAT (grid+TSclim,Rmean,TS at boundaries)
    open(50,FILE='../data/INITnew.DAT',FORM='UNFORMATTED')
    read(50) Z,ZZ,DZ,DZZ,ALON,ALAT,DX,DY,H,COR,ANG, &
    ART,ARU,ARV,FSM,DUM,DVM,T,S,TCLIM,SCLIM, &
    RRMEAN,TBWO,SBWO,TBEO,SBEO,TBSO,SBSO,TBNO,SBNO
    close(50)

    open(32,file='bac3dannMed10fromMed20.bin',form='unformatted')
    do k=1,kb-1
        read(32)(bac3d(i,1,k),i=1,im*jm)
    enddo
    close(32)

!---open NetCDF files
    CFILE='bac3dAnnFromMed20.nc'
    istatus=NF_CREATE(CFILE(1:(INDEX(CFILE,' '))-1),NF_SHARE,NCID1)!create dataset
    write(*,*)'OK',NCID1
    istatus = NF_DEF_DIM(NCID1, 'z',kb-1, KLEVID)
    istatus = NF_DEF_DIM(NCID1, 'y', JM, JMID)
    istatus = NF_DEF_DIM(NCID1, 'x', IM, IMID)
    write(*,*)'OK',KLEVID,JMID,IMID
    call write2dnetcdf(im,jm,kb-1, &
    ncid1,IMID,JMID,KLEVID,bac3d,4,'bact')
    istatus = NF_CLOSE (NCID1)
    open(12,file='bac3dannZfromMed20.bin',form='unformatted')
    call SIGMA2CART(bac3dz, bac3dav, -1.e+10)
    do k=1,klev
        write(12)(bac3dz(i,1,k),i=1,im*jm)
    enddo
    close(12)
    stop
end program





!----------------------------------------------------------------------
subroutine write3dnetcdf(im,jm,kb,nstates,ncid,imid,jmid,klevid,stateid,var,ldf1,varname1)
!! Write 3D data in netCDF file.
    include "netcdf.inc"
    integer :: JMID,IMID,KLEVID,STATEID,VARID
    integer :: VAR3DDIMS(4)
    real :: var(im,jm,kb-1,nstates)
    character varname1*50

    istatus = NF_REDEF(NCID)
    call define_variable_shape3D(VAR3DDIMS,IMID,JMID,KLEVID,STATEID)
    write(*,*)ncid,imid,jmid,klevid,STATEID
    istatus = NF_DEF_VAR(NCID,varname1(1:ldf1),NF_FLOAT, 4, &
    VAR3DDIMS, VARID)
    istatus = NF_ENDDEF(ncid, iret)
    istatus = nf_put_var_real(ncid,VARID,var)
    call check_err(istatus)
!    istatus = NF_CLOSE (NCID)
    write(*,*)'3-d Netcdf'
    return
end subroutine write3dnetcdf




subroutine write2dnetcdf(im,jm,nstates,ncid,IMID,JMID,STATEID,var,ldf1,varname1)
!! Write 2D data in netCDF file.
    include "netcdf.inc"
    integer :: JMID,IMID,VARID,STATEID
    integer :: VAR2DDIMS(3)
    real :: var(im,jm,nstates)
    character varname1*50

    istatus = NF_REDEF(NCID)
    call define_variable_shape2D(VAR2DDIMS,IMID,JMID,STATEID)
    istatus = NF_DEF_VAR(NCID,varname1(1:ldf1),NF_FLOAT,3,VAR2DDIMS,VARID)
    istatus = NF_ENDDEF(ncid, iret)
    istatus = nf_put_var_real(ncid,VARID,var)
    call check_err(istatus)
!    istatus = NF_CLOSE (NCID)
    write(*,*)'2-d Netcdf'
    return
end subroutine write2dnetcdf
!----------------------------------------------------------------------




!----------------------------------------------------------------------
subroutine define_variable_shape4D(IVALDIMS,IMID,JMID,KLEVID,NTIMESID)
!! ??
    integer :: IVALDIMS(4)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    IVALDIMS(3) = KLEVID
    IVALDIMS(4) = NTIMESID
    return
end subroutine define_variable_shape4D


subroutine define_variable_shape3D(IVALDIMS,IMID,JMID,KLEVID,ISTATEID)
!! ??
    integer :: IVALDIMS(4)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    IVALDIMS(3) = KLEVID
    IVALDIMS(4) = ISTATEID
    return
end subroutine define_variable_shape3D



subroutine define_variable_shape2D(IVALDIMS,IMID,JMID,ISTATEID)
!! ??
    integer :: IVALDIMS(3)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    IVALDIMS(3) = ISTATEID
    return
end subroutine define_variable_shape2D



subroutine define_variable_shape1D(IVALDIMS,KLEVID)
!! ??
    integer :: IVALDIMS(1)
    IVALDIMS(1) = KLEVID
    return
end subroutine define_variable_shape1D
!----------------------------------------------------------------------




subroutine define_variable_shape2D(IVALDIMS,IMID,JMID,ISTATEID)
!! ??
    integer :: IVALDIMS(3)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    IVALDIMS(3) = ISTATEID
    return
end subroutine define_variable_shape2D





subroutine define_variable_shape1D(IVALDIMS,KLEVID)
!! ??
    integer :: IVALDIMS(1)
    IVALDIMS(1) = KLEVID
    return
end subroutine define_variable_shape1D


!----------------------------------------------------------------------
subroutine FINDTIME(TIME,DATE)
!! ??
    character DATE*12

    read(DATE,5) IYEAR,MONTH,IDAY,IHOUR,MINUT
    5 format(I4.2,4I2.2)
    JDN = IDAY - 32075 &
    + 1461 * (IYEAR + 4800 + (MONTH - 14) / 12) / 4 &
    + 367 * (MONTH - 2 - (MONTH -14) / 12 * 12) / 12 &
    - 3 * ((IYEAR + 4900 + (MONTH - 14) / 12) / 100) / 4
    TIME = FLOAT(JDN) + FLOAT(IHOUR)/24.
    return
end subroutine FINDTIME





subroutine check_err(iret)
!! ??
    integer :: iret
    include 'netcdf.inc'
    if (iret /= NF_NOERR) then
        print *, nf_strerror(iret)
        stop
    endif
end subroutine check_err




!---Initialization of the 3d<->1d mappings----------------
subroutine setup_3to1_mapping(fsm, im, jm, kb)
!! Initialization of the 3d<->1d mappings
    implicit none
    include 'sesame_params.h'
    include 'mappings.h'
    integer :: i,j,k,iw
    integer :: im,jm,kb
    real :: fsm(im,jm)
    data init_map/0/

!---set init_map to 1 (the mapXtoX functions check this for error checking)
    init_map = 1
    iw = 0
    do k = 1, kb-1
        do j = 1, jm
            do i = 1, im
                if( fsm(i,j) > 0. ) then
                    iw = iw + 1
                    i3d(iw) = i
                    j3d(iw) = j
                    k3d(iw) = k
                    if(i == 320 .AND. j == 37 .AND. k == 1) &
                    write(*,*)'(320,37)=',iw
                endif
            enddo
        enddo
    enddo
    if ( iw /= N_COMPg ) then
        write(*,*) 'inconsistency in water points for ERSEM/POM'
        write(*,*) 'iw = ', iw, ' n_comp = ', N_COMPg
        stop
    endif
    write(*,*)'MAPPING',iw,N_COMPg
    return
end subroutine setup_3to1_mapping



!---Maps 1d arrays to 3d-------------------
subroutine map1to3(ar1d, ar3d)
!! Maps 1D arrays to 3D
    include 'sesame_params.h'
    include 'grid.h'
    include 'mappings.h'
    real :: ar3d(im,jm,kb), ar1d(N_COMPg)
    integer :: i,j,k,n

    if(init_map == 0) then
        write(*,*) 'error in mappings'
        call exit(1)
    endif
    do n = 1, N_COMPg
        i = i3d(n)
        j = j3d(n)
        k = k3d(n)
        ar3d(i, j, k) = ar1d(n)
    enddo
    return
end subroutine map1to3





subroutine map1to2(ar1d,ar2d)
!! Maps 1D arrays to 2D
    include 'sesame_params.h'
    include 'grid.h'
    include 'mappings.h'
    integer :: i,j,n
    real :: ar2d(im,jm), ar1d(N_UPPER$g)

    do n = 1, N_UPPER$g
        i = i3d(n)
        j = j3d(n)
        ar2d(i, j)=ar1d(n)
    enddo
    return
end subroutine map1to2






subroutine PRXYZ(LABEL,TIME,A,IM,ISKP,JM,JSKP,KB,SCALA)
!! This subroutine writes layers of a 3-D field 
!! TIME=time in days
!! A = array(IM,JM,KB) to be printed
!! ISKP=print skip for I
!! JSKP=print skip for J
!! SCALE=divisor for values of A
    dimension A(IM,JM,KB),NUM(430),pLINE(430),KP(4)
    character LABEL*(*)
    data KE,KP /4, 1,2,24,25 /
    data ZERO /1.E-10/

    SCALE=SCALA
    if (SCALE > ZERO) goto 150
    AMX=ZERO
    do 140  KM=1,KE
        K=KP(KM)
        do 140 J=1,JM,JSKP
            do 140 I=1,IM,ISKP
                AMX=MAX(ABS(A(I,J,K)),AMX)
    140 enddo
    if(AMX == 0.) THEN
        SCALEI=0.
        goto 165
    endif
    SCALE=10.E0**(INT(LOG10(AMX)+1.E2)-103)
    150 continue
    SCALEI=1.E0/SCALE
    165 continue
    write(6,160) LABEL
    160 format(1X,A40)
    write(6,170) TIME,SCALE
    170 format(' TIME = ',F9.4,' DAYS     MULTIPLY ALL VALUES BY',1PE10.3)
    do 180 I=1,IM
        NUM(I)=I
    180 enddo

    do 500 KM=1,KE
        K=KP(KM)
        write(6,190) K
        190 format(3X,/7H LAYER ,I2)
        IB=1
        200 continue
        IE=IB+23*ISKP
        if(IE > IM) IE=IM
        write(6,220) (NUM(I),I=IB,IE,ISKP)
        220 format(/,2X,24I5,/)
        do 260 J=1,JM,JSKP
            JWR=JM+1-J
            do 230 I=IB,IE,ISKP
                pLINE(I)=SCALEI*A(I,JWR,K)
            230 enddo
            write(6,240) JWR,(pLINE(I),I=IB,IE,ISKP)
            240 format(1X,I3,24(F5.1))
        260 enddo
        write(6,280)
        280 format(//)
        if(IE >= IM) goto 500
        IB=IB+24*ISKP
        goto 200
    500 enddo
    return
end subroutine PRXYZ

