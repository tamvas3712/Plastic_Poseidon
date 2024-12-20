module poseidon_common_subroutines
!write "use poseidon_common_subroutines" where you need these subroutines
!Implicit NONE

contains

!----------------find time and check errors---------------------


subroutine findtime0(time0,year,iihour,iiday,iimo,iiyr)
!! find time
    dimension month_o(12),month_lp(12),month(0:12)
    data month_o/31,59,90,120,151,181,212,243,273,304,334,365/
    data month_lp/31,60,91,121,152,182,213,244,274,305,335,366/
    year=float(iiyr)
    year_days = 365.
    if (mod(year,4.) == 0.) year_days = 366.
    month(0) = 0
    do k = 1,12
        month(k) = month_o(k)
    end do
    if (year_days == 366.) then
        do k = 1,12
            month(k) = month_lp(k)
        end do
    end if
    time0=float(month(iimo-1)+iiday) + float(iihour)/24.
    return
end subroutine findtime0



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


!------------------PRint XYs----------------


subroutine PRXY(LABEL,TIME,A,IM,ISKP,JM,JSKP,SCALA)
!! This subroutine writes a 2-D field 
!! TIME=time in days
!! A = array(IM,JM) to be printed
!! ISKP=print skip for I
!! JSKP=print skip for J
!! SCALE=divisor for values of A

!---Implicit half precision (A-H,O-Z)
    dimension A(IM,JM),NUM(430),aLINE(430)
    character LABEL*(*)
    data ZERO /1.E-12/

    SCALE=SCALA
    if (SCALE > ZERO) goto 160
    AMX=ZERO
    do 150 J=1,JM,JSKP
        do 150 I=1,IM,ISKP
            AMX=MAX(ABS(A(I,J)),AMX)
    150 enddo
    if(AMX == 0.) THEN
        SCALEI=0.
        goto 165
    endif
    SCALE=10.E0**(INT(LOG10(AMX)+1.E2)-103)
    160 continue
    SCALEI=1.E0/SCALE
    165 continue
    write(6,170) LABEL
    170 format(1X,A40)
    write(6,180) TIME,SCALE
    180 format(' TIME =',F9.4,' DAYS     MULTIPLY ALL VALUES BY',1PE10.3)
    do 190 I=1,IM
        NUM(I)=I
    190 enddo
    IB=1

    200 continue
    IE=IB+23*ISKP
    if(IE > IM) IE=IM
    write(6,210) (NUM(I),I=IB,IE,ISKP)
    210 format(/,2X,24I5,/)
    do 260 J=1,JM,JSKP
        JWR=JM+1-J
        do 220 I=IB,IE,ISKP
            aLINE(I)=INT(SCALEI*A(I,JWR))
        220 enddo
        write(6,240) JWR,(aLINE(I),I=IB,IE,ISKP)

        240 format(1X,I3,24f5.1)
    260 enddo
    write(6,280)
    280 format(//)
    if(IE >= IM) return
    IB=IB+24*ISKP
    goto 200
end subroutine PRXY






subroutine PRXYZ(LABEL,TIME,A,IM,ISKP,JM,JSKP,KB,SCALA)
!! PRINT in X Y Z coordinates. This subroutine writes layers of a 3-D field 
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



!----------randoms-------------------------


integer function randomf(n1,n2)
    implicit none
    integer :: n1,n2
    real*4 :: r
    call random_number(r)
    randomf=n1+int((n2-n1+1)*r)
end function randomf


integer function randomf1(n1,n2)
    implicit none
    integer :: n1,n2
    real*4 :: r
    call random_number(r)
    randomf1=n1+int((n2-n1+1)*r)
end function randomf1


!-------define variable shapes--------------


subroutine define_variable_shape4D(IVALDIMS,IMID,JMID,KLEVID,NTIMESID)
!! ??
    integer :: IVALDIMS(4)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    IVALDIMS(3) = KLEVID
    IVALDIMS(4) = NTIMESID
    return
end subroutine define_variable_shape4D



subroutine define_variable_shape3D(IVALDIMS,IMID,JMID,KLEVID)
!! ??
    integer :: IVALDIMS(3)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    IVALDIMS(3) = KLEVID
    return
end subroutine define_variable_shape3D



subroutine define_variable_shape2D(IVALDIMS,IMID,JMID)
!! ??
    integer :: IVALDIMS(2)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    return
end subroutine define_variable_shape2D



subroutine define_variable_shape2Dt(IVALDIMS,IMID,JMID,NTIMESID)
!! ??
    integer :: IVALDIMS(3)
    IVALDIMS(1) = IMID
    IVALDIMS(2) = JMID
    IVALDIMS(3) = NTIMESID
    return
end subroutine define_variable_shape2Dt



subroutine define_variable_shape1D(IVALDIMS,KLEVID)
!! ??
    integer :: IVALDIMS(1)
    IVALDIMS(1) = KLEVID
    return
end subroutine define_variable_shape1D



!-------------write xd netcdf--------------------



subroutine write1dnetcdf(kb,ncid,KLEVID,var,ldf1,ldf2,ldf3, &
    varname1,varname2,varunits)
!! Write 1D data in netCDF file.
    include "netcdf.inc"
    integer :: KLEVID,VARID
    integer :: VAR1DDIMS(1)
    real :: var(kb)
    character varname1*50,varname2*50,varunits*50

    istatus = NF_REDEF(NCID)
    call define_variable_shape1D(VAR1DDIMS,KLEVID)
!    write(*,*)varname1(1:ldf1-1),var(1),var(kb),kb
    istatus = NF_DEF_VAR(NCID,varname1(1:ldf1-1),NF_FLOAT,1,VAR1DDIMS,VARID)
    write(*,*)VARID,KLEVID,VAR1DDIMS(1)
    istatus = NF_PUT_ATT_TEXT(NCID,VARID,'long_name',ldf2-1,varname2(1:ldf2-1))
    istatus = NF_PUT_ATT_TEXT(NCID,VARID,'units',ldf3-1,varunits(1:ldf3-1))
    istatus = NF_ENDDEF(ncid, iret)
    istatus = nf_put_var_real(ncid,VARID,var)
    call check_err(istatus)
!    istatus = NF_CLOSE (NCID)
    return
end subroutine write1dnetcdf



subroutine write2dnetcdf(im,jm,ncid,IMID,JMID,var,ldf1,ldf2,ldf3, &
    varname1,varname2,varunits)
!! Write 2D data in netCDF file.
    include "netcdf.inc"
    integer :: JMID,IMID,VARID
    integer :: VAR2DDIMS(2),LAT(1),LON(1)
    real :: var(im,jm)
    character varname1*50,varname2*50,varunits*50

    istatus = NF_REDEF(NCID)
    write(*,*)IMID,JMID
    call define_variable_shape2D(VAR2DDIMS,IMID,JMID)
    istatus = NF_DEF_VAR(NCID,varname1(1:ldf1-1),NF_FLOAT,2,VAR2DDIMS,VARID)
    istatus = NF_PUT_ATT_TEXT(NCID, VARID,'long_name',ldf2-1,varname2(1:ldf2-1))
    istatus = NF_PUT_ATT_TEXT(NCID,VARID,'units',ldf3-1,varunits(1:ldf3-1))
    write(*,*)VARID,IMID,JMID
    write(*,*)varname1(1:ldf1-1)
    write(*,*)varname2(1:ldf2)
    istatus = NF_ENDDEF(ncid, iret)
    istatus = nf_put_var_real(ncid,VARID,var)
    call check_err(istatus)
!    istatus = NF_CLOSE (NCID)
    write(*,*)'2-d Netcdf'
    return
end subroutine write2dnetcdf



subroutine write3dnetcdf(im,jm,kb,ncid,imid,jmid,klevid,var, &
    ldf1,ldf2,ldf3, &
    varname1,varname2,varunits)
!! Write 3D data in netCDF file.
    include "netcdf.inc"
    integer :: JMID,IMID,KLEVID,VARID
    integer :: VAR3DDIMS(3)
    real :: var(im,jm,kb)
    character varname1*50,varname2*50,varunits*50

    istatus = NF_REDEF(NCID)
    call define_variable_shape3D(VAR3DDIMS,IMID,JMID,KLEVID)
    istatus = NF_DEF_VAR(NCID,varname1(1:ldf1-1),NF_FLOAT,2,VAR3DDIMS,VARID)
    istatus = NF_PUT_ATT_TEXT(NCID, VARID,'units',ldf3-1,varunits(1:ldf2-1))
    write(*,*)VARID,IMID,JMID
    write(*,*)varname1(1:ldf1-1)
    write(*,*)varname2(1:ldf2)
    istatus = NF_ENDDEF(ncid, iret)
    istatus = nf_put_var_real(ncid,VARID,var)
    call check_err(istatus)
    istatus = NF_CLOSE (NCID)
    write(*,*)'3-d Netcdf'
    return
end subroutine write3dnetcdf



subroutine write4dnetcdf(im,jm,kb,ncid,imid,jmid,klevid,NTIMESID, &
    var,ldf1,ldf2,kld3, &
    varname1,varname2,varunits)
!! Write 4D data in netCDF file.
    include "netcdf.inc"
    integer :: JMID,IMID,KLEVID,NTIMESID,VARID
    integer :: VAR4DDIMS(4)
    real :: var(im,jm,kb,4)
    character varname1*50,varname2*50,varunits*50

    istatus = NF_REDEF(NCID)
    call define_variable_shape4D(VAR4DDIMS,IMID,JMID,KLEVID,NTIMESID)
    istatus = NF_DEF_VAR(NCID,varname1(1:ldf1-1),NF_FLOAT,4,VAR4DDIMS,VARID)
    istatus = NF_PUT_ATT_TEXT(NCID,VARID,'long_name',ldf2-1,varname2(1:ldf2-1))
    istatus = NF_PUT_ATT_TEXT(NCID, VARID,'units',ldf3-1,varunits(1:ldf2-1))
    istatus = NF_ENDDEF(ncid, iret)
    istatus = nf_put_var_real(ncid,VARID,var)
    call check_err(istatus)
    return
end subroutine write4dnetcdf




end module poseidon_common_subroutines
