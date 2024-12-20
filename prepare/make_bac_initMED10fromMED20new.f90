PROGRAM MAIN
!! Prepare MED10 state variable from EMED restart and MEDATLAS climatology for nutrients

!--COARSE-FINE grids in grid.h
    include 'grid.h'

!--COARSE
    dimension bac3d_coa(im2,jm2,kb2)
    parameter (xleft=-7., ybot=30.25, xstep=0.05, ystep=0.05)

!--FINE
    dimension bac3d(im,jm,kb)

!--read EMED model grid
    open(550,FILE='../data/INITmed20.DAT',FORM='UNFORMATTED')
    read(550) Z2,ZZ2,DZ2,DZZ2,ALON2,ALAT2,DX2,DY2,H2
    close(550)
    do j=1,jm2
        do i=1,im2
            fsm2(i,j)=1.
            if(h2(i,j) <= 1.)fsm2(i,j)=0.
        enddo
    enddo
!    CALL PRXY(' H-MED10 ',TIME,H2,IM2,1,JM2,1,1.)
!    CALL PRXY(' ALON-MED ',TIME,ALON2,IM2,1,JM2,1,1.E-1)
!    write(*,*)'TEST',(alon2(i,1),i=1,10)
!    write(*,*)'TEST',(alat2(1,j),j=1,10)

!--read MED10 model grid
    open(550,FILE='../data/INITmed10.DAT',FORM='UNFORMATTED')
    read(550) Z,ZZ,DZ,DZZ,ALON,ALAT,DX,DY,H
    close(550)
    CALL PRXY(' TOPO',TIME,H,IM,1,JM,1,1.)
!    CALL PRXY(' FSM-NAEG ',TIME,FSM,IM,1,JM,1,1.)
    do j=1,jm
        do i=1,im
            fsm(i,j)=1.
            if(h(i,j) <= 1.)fsm(i,j)=0.
        enddo
    enddo

!--read MED10 Bacteria
    open(11,file='../data/bac3dannMed20.bin',form='unformatted')
    do k=1,kb2-1
        read(11)(bac3d_coa(i,1,k),i=1,im2*jm2)
    enddo
    close(11)
    CALL PRXYZ('VAR-COARSE',TIME,bac3d_coa,IM2,1,JM2,1,KB2,1.E0)
     
!--interpolate into MED10 grid
    call MAPDATA(bac3d_coa,bac3d,xleft,ybot,xstep,ystep)
    CALL PRXYZ('VAR',TIME,bac3d,IM,1,JM,1,KB,1.E0)

!--write MED20  bacteria
    open(11,file='../data/bac3dannMed10.bin',form='unformatted')
    do k=1,kb-1
        write(11)(bac3d(i,1,k),i=1,im*jm)
    enddo
    close(11)
end PROGRAM

     







subroutine MAPDATA(a,aa,xleft,ybot,xstep,ystep)
!! Write map data ??          
    include 'grid.h'
    dimension A(IM2,JM2,KB2),AA(IM,JM,KB)
    dimension tpr(kb2),zm2(kb2)
    xland=999.

    write(*,*)'OK'


    do 10 j=1,jm
        do 10 i=1,im
!            i=102
!            j=307
            do k=1,kb
                aa(i,j,k)=xland
            enddo
            if(fsm(i,j) == 0.) go to 10
            x=alon(i,j)
            y=alat(i,j)

        !----locate grid point
            ii=(x-xleft)/xstep+1.
            jj=(y-ybot)/ystep+1.

        !---- find surrounding points--------------------
            x1=xleft+(ii-1)*xstep
            y1=ybot+(jj-1)*ystep
            x2=x1
            y2=y1+ystep
            x3=x1+xstep
            y3=y2
            x4=x3
            y4=y1
!            write(*,*)'MAP',i,j,h(i,j)
            do 15 k=1,kb-2
                imean=0
                fmean=0.
                zm=(-1)*zz(k)*h(i,j)*fsm(i,j)

            !--interpolation for point 1
                do kk = 1, kb2-2
                    tpr(kk)=a(ii,jj,kk)
                enddo
!                write(*,*)'POINT 1',ii,jj
                do kk=1,kb2-2
                    zm2(kk)=(-1)*zz2(kk)*h2(ii,jj)*fsm2(ii,jj)
            !        write(*,*)zm2(kk),tpr(kk),zz2(kk),h2(ii,jj),fsm2(ii,jj)
            !        write(*,*)zm2(kk),tpr(kk),h2(ii,jj)
                enddo
                CALL SINTER(zm2,tpr,zm,ti,kb2-1,1,xland)
                f1=ti
                if(f1 /= xland)then
                    fmean=fmean+f1
                    imean=imean+1
                endif
                ti=0.
        !        write(*,*)f1

            !--interpolation for point 2
                do kk = 1, kb2-2
                    tpr(kk)=a(ii,jj+1,kk)
                enddo
        !        write(*,*)'POINT 2',ii,jj+1
                do kk=1,kb2-2
                    zm2(kk)=(-1)*zz2(kk)*h2(ii,jj+1)*fsm2(ii,jj+1)
            !        write(*,*)zm2(kk),tpr(kk)
                enddo
                CALL SINTER(zm2,tpr,zm,ti,kb2-1,1,xland)
                f2=ti
                if(f2 /= xland)then
                    fmean=fmean+f2
                    imean=imean+1
                endif
                ti=0.
        !        write(*,*)f2

            !--interpolation for point 3
                do kk = 1, kb2-2
                    tpr(kk)=a(ii+1,jj+1,kk)
                enddo
        !        write(*,*)'POINT 3',ii,jj+1
                do kk=1,kb2-2
                    zm2(kk)=(-1)*zz2(kk)*h2(ii+1,jj+1)* &
                    fsm2(ii+1,jj+1)
                !       write(*,*)zm2(kk),tpr(kk)
                enddo

                CALL SINTER(zm2,tpr,zm,ti,kb2-1,1,xland)
                f3=ti
                if(f3 /= xland)then
                    fmean=fmean+f3
                    imean=imean+1
                endif
                ti=0.
            !      write(*,*)f3
            !k--interpolation for point 4

                do kk = 1, kb2-2
                    tpr(kk)=a(ii+1,jj,kk)
                enddo
        !        write(*,*)'POINT 4',ii+1,jj
                do kk=1,kb2-2
                    zm2(kk)=(-1)*zz2(kk)*h2(ii+1,jj)*fsm2(ii+1,jj)
            !        write(*,*)zm2(kk),tpr(kk)
                enddo
                CALL SINTER(zm2,tpr,zm,ti,kb2-1,1,xland)
                f4=ti
                if(f4 /= xland)then
                    fmean=fmean+f4
                    imean=imean+1
                endif
                ti=0.
        !        write(*,*)f4
                fmean=fmean/float(imean)
                ff=amin1(f1,f2,f3,f4)
                if(ff == xland)goto 15
                if(f1 == xland)f1=fmean
                if(f2 == xland)f2=fmean
                if(f3 == xland)f3=fmean
                if(f4 == xland)f4=fmean
                CALL BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
                aa(i,j,k)=f
        !        write(*,*)x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f
            15 enddo

        !--now the benthic variables
            imean=0
            fmean=0.
            if(fsm2(ii,jj) == 1.)then
                f1=a(ii,jj,kb2-1)
                fmean=fmean+f1
                imean=imean+1
            else
                f1=xland
            endif
            if(fsm2(ii,jj+1) == 1.)then
                f2=a(ii,jj+1,kb2-1)
                fmean=fmean+f2
                imean=imean+1
            else
                f2=xland
            endif
            if(fsm2(ii+1,jj+1) == 1.)then
                f3=a(ii+1,jj+1,kb2-1)
                fmean=fmean+f3
                imean=imean+1
            else
                f3=xland
            endif
            if(fsm2(ii+1,jj) == 1.)then
                f4=a(ii+1,jj,kb2-1)
                fmean=fmean+f4
                imean=imean+1
            else
                f4=xland
            endif

            fmean=fmean/float(imean)

            ff=amin1(f1,f2,f3,f4)
            if(ff == xland)goto 10
            if(f1 == xland)f1=fmean
            if(f2 == xland)f2=fmean
            if(f3 == xland)f3=fmean
            if(f4 == xland)f4=fmean
            CALL BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
            aa(i,j,kb-1)=f
        10 enddo
     
!---check for correct interpolation
    do i = 1,im
        do j = 1,jm
            do k=1,kb-2
                if (fsm(i,j) == 1.) then
                    if (aa(i,j,k) == xland) then
                !        print *,'ERROR IN INTERP',I,J,-zz(k)*h(i,j)
                !        print *, 'I,J FOR ERROR ',i,j
                        call correct(aa,i,j,k,xland)
                    end if
                    if (aa(i,j,k) == xland) then
                !        print *,'ERROR IN SPATIAL INTERP',I,J,K,fsm(i,j)
                !        aa(i,j,k)=999.
                    endif
                else
                    aa(i,j,k)=0.
                end if
            end do
        end do
    enddo

    do i = 1,im
        do j = 1,jm
            if(i == 1)write(*,*)'TESTGIB0',j,aa(i,j,kb-1)
            if (fsm(i,j) == 1.) then
                if (aa(i,j,kb-1) == xland) then
                    if(i == 1)write(*,*)'TESTGIB1',j,aa(i,j,kb-1)
                    call correct_ben(aa,i,j,kb-1,xland)
                    if(i == 1)write(*,*)'TESTGIB2',j,aa(i,j,kb-1)
                endif
                if (aa(i,j,kb-1) == xland) then
            !        print *,'ERROR IN SPATIAL INTERP-BENTH',I,J,K,fsm(i,j)
            !        aa(i,j,k)=999.
                endif
            else
                aa(i,j,kb-1)=0.
            endif
        end do
    end do

    return
end subroutine MAPDATA








subroutine correct(aux,ie,je,ke,xland)
!! Perform error correction ?? 
    include 'grid.h'

    dimension aux(im,jm,kb),tpr(kb),zm(kb)
    zm1=(-1)*zz(ke)*h(ie,je)*fsm(ie,je)
    lr = 1

    50 continue

    ist = ie -lr
    iend = ie + lr
    jst = je - lr
    jend = je + lr
    if (ist < 1) ist = 1
    if (iend > im) iend = im
    if (jst < 1) jst = 1
    if (jend > jm) jend = jm

    cc = 0.0
    xmea = 0.0
    do 100 i = ist,iend
        do 100 j = jst,jend
            do k = 1, kb-2
                tpr(k)=aux(i,j,k)
            enddo
            do k=1,kb-2
                zm(k)=(-1)*zz(k)*h(i,j)*fsm(i,j)
            enddo
            CALL SINTER(zm,tpr,zm1,ti,kb-1,2,xland)
            f1=ti
            ti=0.
            if (f1 /= xland) then
                cc = cc + 1.
                xmea = xmea + f1
            endif
    100 enddo

    if (cc /= 0.) then
        aux(ie,je,ke) = xmea/cc
    else
        lr = lr + 1
        if(lr < 20)goto 50
    endif

!    print *,'ERROR CORRECTION FOR ',IE,JE,'WITH LR=',lr

    return
end subroutine correct






subroutine correct_ben(aux,ie,je,ke,xland)
!! Correct benthic variables ??
    include 'grid.h'
    dimension aux(im,jm,kb)

    lr = 1

    50 continue

    ist = ie -lr
    iend = ie + lr
    jst = je - lr
    jend = je + lr
    if (ist < 1) ist = 1
    if (iend > im) iend = im
    if (jst < 1) jst = 1
    if (jend > jm) jend = jm

    cc = 0.0
    xmea = 0.0
    do 100 i = ist,iend
        do 100 j = jst,jend
            if (aux(i,j,ke) /= xland .AND. fsm(i,j) == 1.) then
                cc = cc + 1.
                xmea = xmea + aux(i,j,ke)
                if(ie == 1)write(*,*)j,xmea,cc
            end if
    100 enddo

    if (cc /= 0.) then
        aux(ie,je,ke) = xmea/cc
    else
        lr = lr + 1
        if(lr < 20)goto 50
    end if

!       print *,'ERROR CORRECTION FOR ',IE,JE,'WITH LR=',lr

    return
end subroutine correct_ben






subroutine correct2d(aux,ie,je,xland)
!! Correct 2D variables ??
    include 'grid.h'
    dimension aux(im,jm)

    lr = 1

    50 continue

    ist = ie -lr
    iend = ie + lr
    jst = je - lr
    jend = je + lr
    if (ist < 1) ist = 1
    if (iend > im) iend = im
    if (jst < 1) jst = 1
    if (jend > jm) jend = jm

    cc = 0.0
    xmea = 0.0
    do 100 i = ist,iend
        do 100 j = jst,jend
            if (aux(i,j) /= xland .AND. fsm(i,j) == 1.) then
                cc = cc + 1.
                xmea = xmea + aux(i,j)
            end if
    100 enddo

    if (cc /= 0.) then
        aux(ie,je) = xmea/cc
    else
        lr = lr + 1
        if(lr < 50)goto 50
    end if

!       print *,'ERROR CORRECTION FOR ',IE,JE,'WITH LR=',lr

    return
end subroutine correct2d






subroutine BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
!!  Bilinear interpolation subroutine.
!!  (Xi,Yi,fi) = data grid & values surounding model point (x,y)
!!  f = interpolated value at the model grid point.

    a1=x1-x2+x3-x4
    a2=-x1+x4
    a3=-x1+x2
    a4=x1-x
    b1=y1-y2+y3-y4
    b2=-y1+y4
    b3=-y1+y2
    b4=y1-y
    A=a3*b1-a1*b3
    B=b2*a3+b1*a4-a1*b4-a2*b3
    C=-a2*b4+a4*b2
    if(ABS(A*C) > 0.002*B**2) then
        t=(-B-sqrt(B*B-4.*A*C))/(2.*A)
    else
        t=C/ABS(B)
    endif
    10 continue
    A=a2*b1-a1*b2
    B=b3*a2+b1*a4-a1*b4-a3*b2
    C=-a3*b4+a4*b3
    if(ABS(A*C) > 0.002*B**2) then
        s=(-B+sqrt(B*B-4.*A*C))/(2.*A)
    else
        s=-C/ABS(B)
    endif
    20 continue
    f=f1*(1.-t)*(1.-s)+f2*t*(1.-s)+f3*s*t+f4*(1.-t)*s
    return
end subroutine BLINT








subroutine SINTER(X,A,Y,B,M,imode,xland)
!! This routine linearly interpolates and extrapolates an array B, with M being the number of points in X and A.
    iflag=0

!-- no extrapolation--
    if(Y < X(1)) then
        B=A(1)
        iflag=1
    !         write(*,*)'TOP'
    endif
!    if(Y.GT.X(1)) B=A(1)+((A(1)-A(2))/(X(1)-X(2)))*(Y-X(1))

    go to (16,10,15) imode

    10 continue
    if(Y > X(M-1) .AND. X(M-1) > 0.) then
        B=A(M-1)
!        B=A(M-2) -(A(M-2)-A(M-1))*(X(M-2)-Y)/(X(M-2)-X(M-1))
        iflag=1
!        write(*,*)'BOT'
    endif
    goto 16

    15 continue
    if(Y > X(M-1) .AND. X(M-1) < 2000. .AND. X(M-1) > 0.) then
        B=A(M-1)
        iflag=1
    endif
    goto 16

    16 continue

!---interpolation cases
    NM=M-2
    do 20 J=1,NM
        if(A(J) == xland .OR. A(J+1) == xland)goto 20
        if(Y >= X(J) .AND. Y <= X(J+1))then
            B=A(J) -(A(J)-A(J+1))*(X(J)-Y)/(X(J)-X(J+1))
            iflag=1
    !        write(*,*)'BETWEEN'
        endif
    20 enddo
    if(iflag == 0)then
        b=xland
!        write(*,*)'NO-POINT'
    endif
    if(imode /= 2)then
        do k=1,M
            A(k)=0.
        enddo
    endif
    return
end subroutine SINTER







subroutine RSINTER(ZM,ZLEV,AIN,AZOUT,KLEV,XLAND)
!! Another type of interpolation ??
!!        Z(KB) must be ascending 
!!        AIN(KB) given function 
!!        AZOUT(KLEV) found by linear interpolation and extrapolation
!!        ZLEV(KLEV) the desired depths 

    parameter(KB=25)
    dimension ZM(KB),ZLEV(KLEV),AIN(KB),AZOUT(KLEV)

!---extrapolation cases
    do 30 k=1,klev
        if(zm(1) > zlev(k))then
            AZOUT(k)=AIN(1)
        endif
        if(zm(kb-1) < zlev(k))then
            AZOUT(k)=0.
        endif
    30 enddo

!---interpolation cases
    do 10 k=1,klev
        do 20 kk=1,kb-2
            if(zm(kk) <= zlev(k) .AND. zm(kk+1) >= zlev(k)) then
                AZOUT(k)=(AIN(kk+1)*(zlev(k)-zm(kk))+ &
                AIN(kk)*(zm(kk+1)-zlev(k)))/(zm(kk+1)-zm(kk))
            endif
        20 enddo
    10 enddo

    return
end subroutine RSINTER







subroutine PRXY(LABEL,TIME,A,IM,ISKP,JM,JSKP,SCALA)
!! This subroutine writes a 2D field ??
!!      TIME= time in days
!!      A = ARRAY(IM,JM) to be printed
!!      ISKP= print skip for I
!!      JSKP= print skip for J
!!      SCALE= divisor for values of A

!---implicit half precision (A-H,O-Z)
    dimension A(IM,JM),NUM(980),pLINE(980)
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
    write(*,170) LABEL
    170 format(1X,A40)
    write(*,180) TIME,SCALE
    180 format(' TIME =',F9.4,' DAYS     MULTIPLY ALL VALUES BY',1PE9.2)
    do 190 I=1,IM
        NUM(I)=I
    190 enddo
    IB=1

    200 continue
    IE=IB+23*ISKP
    if(IE > IM) IE=IM
    write(*,210) (NUM(I),I=IB,IE,ISKP)
    210 format(/,2X,24I5,/)
    do 260 J=1,JM,JSKP
        JWR=JM+1-J
        do 220 I=IB,IE,ISKP
            pLINE(I)=(SCALEI*A(I,JWR))
        220 enddo
        write(*,240) JWR,(pLINE(I),I=IB,IE,ISKP)
        240 format(1X,I3,24(F5.1))
    260 enddo
    write(*,280)
    280 format(//)
    if(IE >= IM) return
    IB=IB+24*ISKP
    goto 200
end subroutine PRXY







subroutine PRXYZ(LABEL,TIME,A1,IM,ISKP,JM,JSKP,KB1,SCALA)
!! This writes horizontal layers of a 3D field 
!     TIME=TIME IN DAYS
!     A = ARRAY(IM,JM,KB) TO BE PRINTED
!     ISKP=PRINT SKIP FOR I
!     JSKP=PRINT SKIP FOR J
!     SCALE=DIVISOR FOR VALUES OF A1

    parameter (KE=5)
    dimension A1(IM,JM,KB1),NUM(950),pLINE(950),KP(KE)
    character LABEL*(*)
    data KP /1,2,3,17,18/
    data ZERO /1.E-10/
    KP(4)=KB1-2
    KP(5)=KB1-1
    SCALE=SCALA
    if (SCALE > ZERO) goto 150
    AMX=ZERO
    do 140  KM=1,KE
        K=KP(KM)
        do 140 J=1,JM,JSKP
            do 140 I=1,IM,ISKP
                AMX=MAX(ABS(A1(I,J,K)),AMX)
    140 enddo
    if(AMX == 0.) THEN
        SCALEI=0.
        goto 165
    endif
    SCALE=10.E0**(INT(LOG10(AMX)+1.E2)-103)
    150 continue
    SCALEI=1.E0/SCALE
    165 continue
    write(*,160) LABEL
    160 format(1X,A40)
    write(*,170) TIME,SCALE
    170 format(' TIME = ',F9.4,' DAYS     MULTIPLY ALL VALUES BY',1PE10.3)
    do 180 I=1,IM
        NUM(I)=I
    180 enddo
    do 500 KM=1,KE
        K=KP(KM)
        write(*,190) K
        190 format(3X,/7H LAYER ,I2)
        IB=1
    
        200 continue
        IE=IB+23*ISKP
        if(IE > IM) IE=IM
        write(*,220) (NUM(I),I=IB,IE,ISKP)
        220 format(/,2X,24I5,/)
        do 260 J=1,JM,JSKP
            JWR=JM+1-J
            do 230 I=IB,IE,ISKP
                pLINE(I)=SCALEI*A1(I,JWR,K)
            230 enddo
            write(*,240) JWR,(pLINE(I),I=IB,IE,ISKP)
            240 format(1X,I3,24(F5.1))
        260 enddo
        write(*,280)
        280 format(//)
        if(IE >= IM) goto 500
        IB=IB+24*ISKP
        goto 200
    500 enddo
    return
end subroutine PRXYZ




