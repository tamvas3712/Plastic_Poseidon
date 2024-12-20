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

subroutine get_time
!!  Return the model time
    implicit none
    include 'pom.h'
    time=dti*float(iint)/86400.e0+time0
    if(iint >= iswtch) iprint=nint(prtd2*24.e0*3600.e0/dti)
    return
end subroutine get_time


SUBROUTINE CLOCK4(time,year,mp,mm,xp,xm)
!! CLOCK4 is for seasonal time interpolation , CLOCK12 is for monthly
    DIMENSION DATES(4),DATESN(4),DATESL(4)
    DATA DATESN/45.,135.5,227.,319./
    DATA DATESL/45.5,136.5,228.,320./
          
    year_days = 365.
    do k = 1,4
        dates(k) = datesn(k)
    end do
    if (mod(year,4.) == 0.) then
        year_days = 366.
        do k = 1,4
            dates(k) = datesl(k)
        end do
    end if

    XTIME=MOD(TIME,year_days)

    if (xtime > dates(4)) then
        mm = 4
        mp = 1
        daysp = dates(1)+year_days
        daysm = dates(4)
        goto 160
    end if

    do 100 k = 1,4
        if (xtime <= dates(k)) then
            mp = k
            mm = k-1
            goto 150
        end if

    100 END DO

    150 continue

    daysp = dates(mp)
    if (mm /= 0) daysm = dates(mm)

    if (mm == 0) then
        mm = 4
        daysm = dates(mm)
        xtime = year_days + xtime
        daysp = daysp + year_days
    end if
          

    160 continue

    dist = daysp - daysm

    xp = 1. - (daysp - xtime)/dist
    xm = 1.- xp

    return
END SUBROUTINE CLOCK4




SUBROUTINE CLOCK4_360(time,mp,mm,xp,xm)
!! clock
    DIMENSION DATES(4)
    DATA DATES/45.,135.,225.,315./
    year_days = 360.
    XTIME=MOD(TIME,year_days)

    if (xtime > dates(4)) then
        mm = 4
        mp = 1
        daysp = dates(1)+year_days
        daysm = dates(4)
        goto 160
    end if

    do 100 k = 1,4
        if (xtime <= dates(k)) then
            mp = k
            mm = k-1
            goto 150
        end if

    100 END DO

    150 continue

    daysp = dates(mp)
    if (mm /= 0) daysm = dates(mm)

    if (mm == 0) then
        mm = 4
        daysm = dates(mm)
        xtime = year_days + xtime
        daysp = daysp + year_days
    end if
          

    160 continue

    dist = daysp - daysm

    xp = 1. - (daysp - xtime)/dist
    xm = 1.- xp
    return
END SUBROUTINE CLOCK4_360




SUBROUTINE PRXY(LABEL,TIME,A,IM,ISKP,JM,JSKP,SCALA)
!!      THIS WRITES A 2-D FIELD
!!      TIME=TIME IN DAYS
!!      A = ARRAY(IM,JM) TO BE PRINTED
!!      ISKP=PRINT SKIP FOR I
!!      JSKP=PRINT SKIP FOR J
!!      SCALE=DIVISOR FOR VALUES OF A

!   IMPLICIT HALF PRECISION (A-H,O-Z)
    DIMENSION A(IM,JM),NUM(430),aLINE(430)
    CHARACTER LABEL*(*)
    DATA ZERO /1.E-12/

    SCALE=SCALA
    IF (SCALE > ZERO) GO TO 160
    AMX=ZERO
    DO 150 J=1,JM,JSKP
        DO 150 I=1,IM,ISKP
            AMX=MAX(ABS(A(I,J)),AMX)
    150 END DO
    IF(AMX == 0.) THEN
        SCALEI=0.
        GOTO 165
    ENDIF
    SCALE=10.E0**(INT(LOG10(AMX)+1.E2)-103)
    160 CONTINUE
    SCALEI=1.E0/SCALE
    165 CONTINUE
    WRITE(6,170) LABEL
    170 FORMAT(1X,A40)
    WRITE(6,180) TIME,SCALE
    180 FORMAT(' TIME =',F9.4,' DAYS     MULTIPLY ALL VALUES BY',1PE10.3)
    DO 190 I=1,IM
        NUM(I)=I
    190 END DO
    IB=1

    200 CONTINUE
    IE=IB+23*ISKP
    IF(IE > IM) IE=IM
    WRITE(6,210) (NUM(I),I=IB,IE,ISKP)
    210 FORMAT(/,2X,24I5,/)
    DO 260 J=1,JM,JSKP
        JWR=JM+1-J
        DO 220 I=IB,IE,ISKP
            aLINE(I)=INT(SCALEI*A(I,JWR))
        220 END DO
        WRITE(6,240) JWR,(aLINE(I),I=IB,IE,ISKP)

        240 FORMAT(1X,I3,24f5.1)
    260 END DO
    WRITE(6,280)
    280 FORMAT(//)
    IF(IE >= IM) RETURN
    IB=IB+24*ISKP
    GO TO 200
END SUBROUTINE PRXY



SUBROUTINE PRXYZ(LABEL,TIME,A,IM,ISKP,JM,JSKP,KB,SCALA)
! >>>

!     THIS WRITES HORIZONTAL LAYERS OF A 3-D FIELD
!     TIME=TIME IN DAYS
!     A = ARRAY(IM,JM,KB) TO BE PRINTED
!     ISKP=PRINT SKIP FOR I
!     JSKP=PRINT SKIP FOR J
!     SCALE=DIVISOR FOR VALUES OF A

    DIMENSION A(IM,JM,KB),NUM(430),pLINE(430),KP(4)
    CHARACTER LABEL*(*)
    DATA KE,KP /4, 1,2,24,25 /
    DATA ZERO /1.E-10/

    SCALE=SCALA
    IF (SCALE > ZERO) GO TO 150
    AMX=ZERO
    DO 140  KM=1,KE
        K=KP(KM)
        DO 140 J=1,JM,JSKP
            DO 140 I=1,IM,ISKP
                AMX=MAX(ABS(A(I,J,K)),AMX)
    140 END DO
    IF(AMX == 0.) THEN
        SCALEI=0.
        GOTO 165
    ENDIF
    SCALE=10.E0**(INT(LOG10(AMX)+1.E2)-103)
    150 CONTINUE
    SCALEI=1.E0/SCALE
    165 CONTINUE
    WRITE(6,160) LABEL
    160 FORMAT(1X,A40)
    WRITE(6,170) TIME,SCALE
    170 FORMAT(' TIME = ',F9.4,' DAYS     MULTIPLY ALL VALUES BY',1PE10.3)
    DO 180 I=1,IM
        NUM(I)=I
    180 END DO

    DO 500 KM=1,KE
        K=KP(KM)
        WRITE(6,190) K
        190 FORMAT(3X,/7H LAYER ,I2)
        IB=1
    
        200 CONTINUE
        IE=IB+23*ISKP
        IF(IE > IM) IE=IM
        WRITE(6,220) (NUM(I),I=IB,IE,ISKP)
        220 FORMAT(/,2X,24I5,/)
        DO 260 J=1,JM,JSKP
            JWR=JM+1-J
            DO 230 I=IB,IE,ISKP
                pLINE(I)=SCALEI*A(I,JWR,K)
            230 END DO
            WRITE(6,240) JWR,(pLINE(I),I=IB,IE,ISKP)
            240 FORMAT(1X,I3,24(F5.1))
        260 END DO
        WRITE(6,280)
        280 FORMAT(//)
        IF(IE >= IM) GO TO 500
        IB=IB+24*ISKP
        GO TO 200
    500 END DO
    RETURN
    END SUBROUTINE PRXYZ

    SUBROUTINE CLOCK12(time,year,mp,mm,xp,xm)

    dimension dates(12),datesn(12),datesl(12)

    data datesn/ 15.5,  45.,   74.5, 105.,  135.5, 166., &
    196.5, 227.5, 258.,  288.5, 319.,  349.5 /

    data datesl/ 15.5,  45.5,  75.5, 106.,  136.5, 167., &
    197.5, 228.5, 259.,  289.5, 320.,  350.5 /
          
!******************************************************************
          
          
    year_days = 365.

    do k = 1,12
        dates(k) = datesn(k)
    end do

    if (mod(yearp,4.) == 0.) then
        year_days = 366.
        do k = 1,12
            dates(k) = datesl(k)
        end do
    end if

    xtime=mod(time,year_days)
    if (xtime > dates(12)) then
        mm = 12
        mp = 1
        daysp = dates(1)+year_days
        daysm = dates(12)
    else
        do k = 1,12
            if (xtime <= dates(k)) then
                mp = k
                mm= k-1
                goto 654
            endif
        enddo
        654 continue
        daysp = dates(mp)
        daysm = dates(mm)
        if (mm == 0) then
            mm = 12
            daysm = dates(mm)
            xtime = year_days + xtime
            daysp = daysp + year_days
        end if
    endif

    dist = daysp - daysm

    xp = 1. - (daysp - xtime)/dist
    xm = 1.- xp
                                                      
          
    return
    END SUBROUTINE CLOCK12

    subroutine findtime0(time0,year,iihour,iiday,iimo,iiyr)
!! Find time in days based on date
! -- Univ. of Athens - Ocean
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


