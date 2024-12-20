program MAIN
!! main program about what ??
    include 'comblk.h'
    include 'fsm3d.h'
!    include 'ecooutput.h'
          
    dimension var3d(im,jm,klev),var2d(im,jm)
    dimension var3dav(im,jm,klev,28),var2dav(im,jm)
    dimension var3dhyav(im,jm,klev,5),var2dbotav(im,jm,5)
    dimension amld(im,jm),amldav(im,jm)
    dimension outmon(4)
    character directory*110,outmon*7,outyr*2,outname*4, &
    outexte*4,fname*110,chara*1
    data outmon/'Jan_Mar','Apr_Jun','Jul_Sep','Oct_Dec'/
    data outname/'erun'/
    data outexte/'.bin'/
    real :: aland
    data aland / -1.0E+10 /
    integer :: ecoidfile1,ecoidfile2
    real :: fsm2(im,jm),h2(im,jm)

!---read MED10 model grid
    open(550,FILE='../data/INITnew.DAT',FORM='UNFORMATTED')
    read(550) Z,ZZ,DZ,DZZ,ALON,ALAT,DX,DY,H
    close(550)

!    open(550,FILE='../data/bath1000Med20.bin',FORM='UNFORMATTED')
!    read(550) (H2(i,1),i=1,im*jm)
!    close(550)

!    do j=1,jm
!        do i=1,im
!            fsm(i,j)=1.
!            fsm2(i,j)=1.
!            if(h(i,j).le.1.)fsm(i,j)=0.
!            if(h2(i,j).ge.500.)fsm2(i,j)=0.
!        enddo
!    enddo

!    open(134,file='../data/topomed20.bin',FORM='UNFORMATTED')
!    write(134)(h(i,1),i=1,im*jm)
!    write(134)(fsm(i,1),i=1,im*jm)
!    close(134)

!    stop

!    open(44,file='../data/chln1popeMed10_2007newmo.bin',form='unformatted')
!    open(44,file='../data/erun2008med20SeasPOS.bin',form='unformatted')
!    open(44,file='../data/erun2010med20monPOSassim.bin',form='unformatted')
    open(44,file='../data/erun2010_2014med20seasPOSoxy23.bin',form='unformatted')
!    open(44,file='../data/zoo2016_2018med20mon2.bin',form='unformatted')

!    open(42,file='../data/chlMed10_2000moDepo.bin',form='unformatted')
!    open(43,file='../data/sst_sss_mld1990_1996medDMIMo.bin',
!    open(42,file='../data/par_euz2017med20monPOStest.bin',form='unformatted')

    do 2 iyear=2010,2014
        directory='../data/'
    !    directory='/home/ktsiaras/BLACKSEA/HindEco/'
        ldf=0
        do ii=1,110
            chara=directory(ii:ii)
            if(chara /= ' ')then
                ldf=ldf+1
            endif
        enddo

        imo=1
        do 5 iseas=1,4
            if(iyear >= 2000)then
                write(outyr,'(i2.2)')iyear-2000
            else
                write(outyr,'(i2.2)')iyear-1900
            endif
            fname=directory(1:ldf)//outname//outmon(iseas)//outyr//outexte
            ldf1=ldf+17
            write(*,*)fname(1:ldf1),outyr
            open(33,file=fname(1:ldf1),form='unformatted',status='old')

            ecoidfile1=33
            ecoidfile2=44
        !    iaver=3*3
            iaver=3
        !    iaver=1
            aver=float(iaver)
        !    nend=9
            nend=3
        !    if(iseas.eq.2)then
        !        do nout=1,3
        !            do nvar=1,20
        !                call read_variable(ecoidfile1, var3d, aland)
        !            enddo
        !            call read_variable_2d(ecoidfile1, var2d, aland)
        !            do nvar=1,5
        !                call read_variable_hydro(ecoidfile1, var3d , aland)
        !            enddo
        !        enddo
        !    endif

        !    if(iseas.eq.2)nend=6
            do 10 nout=1,nend
                do nvar=1,23
                    call read_variable(ecoidfile1, var3d, aland)
                    do k=1,klev
                        do j=1,jm
                            do i=1,im
                                var3dav(i,j,k,nvar)=var3dav(i,j,k,nvar)+var3d(i,j,k)/aver
                            enddo
                        enddo
                    enddo
                enddo
                do nvar=1,1
                    call read_variable_2d(ecoidfile1, var2d, aland)
                    do j=1,jm
                        do i=1,im
                            var2dbotav(i,j,nvar)=var2dbotav(i,j,nvar)+var2d(i,j)/aver
                        enddo
                    enddo
                enddo
                call read_variable_2d(ecoidfile1, var2d, aland)
                write(*,*)'OK'
                do j=1,jm
                    do i=1,im
                        var2dav(i,j)=var2dav(i,j)+var2d(i,j)/aver
                    enddo
                enddo
                do nvar=1,5
                    call read_variable_hydro(ecoidfile1, var3d , aland)
                    if(nvar == 3)t(i,j,k)=var3d(i,j,k)
                    if(nvar == 4)s(i,j,k)=var3d(i,j,k)
                    do k=1,klev
                        do j=1,jm
                            do i=1,im
                                if(nvar == 3)t(i,j,k)=var3d(i,j,k)
                                if(nvar == 4)s(i,j,k)=var3d(i,j,k)
                                var3dhyav(i,j,k,nvar)=var3dhyav(i,j,k,nvar)+var3d(i,j,k)/aver
                            enddo
                        enddo
                    enddo
                    write(*,*)'OK',nvar,var3dhyav(120,100,10,nvar)
                enddo

            !--calculate MLD
                CALL DENS(S,T,RHO,dpt)
                CALL PRXYZ('RHO',TIME,rho,IM,10,JM,10,KB,1.E0)

                do 15 j=1,jm
                    do 15 i=1,im
                        if(fsm(i,j) == 0.)goto 15
                        do k=2,kbm1
                            zm=dpt(k)
                        !   if(abs(t(i,j,1)-t(i,j,k)).ge.0.5)then
                            if(abs(rho(i,j,1)-rho(i,j,k)) >= 0.15)then
                                amld(i,j)=zm
                                goto 15
                            endif
                        enddo
                15 enddo
                do j=1,jm
                    do i=1,im
                        if(amld(i,j) == 0. .AND. fsm(i,j) == 1.)then
                            amld(i,j)=h(i,j)
                        endif
                        if(fsm(i,j) == 0.)then
                            amld(i,j)=aland
                        endif
                    enddo
                enddo
                write(*,*)'MLD',amld(120,100),t(120,100,1),s(120,100,1), &
                rho(120,100,1)
                do j=1,jm
                    do i=1,im
                        amldav(i,j)=amldav(i,j)+amld(i,j)/aver
                    enddo
                enddo

                if(mod(nout,iaver) == 0)then
                    imo=imo+1
                !---calculate also pCO2 at 20 degrees Celsius (write this instead of cco2)
                    do k=1,klev
                        do j=1,jm
                            do i=1,im
                                temp=var3dhyav(i,j,k,1)
                                pco2w=var3dav(i,j,k,27)
                                if(temp > 0.)then
                                    var3dav(i,j,k,28)=pco2w* &
                                !     &     exp(0.0433*(14.9-temp)-4.35*1.E-5*(222.-temp**2.))
                                    exp(0.0433*(13.-temp)-4.35*1.E-5*(169.-temp**2.))
                                !     &     exp(0.0433*(20.-temp)-4.35*1.E-5*(400.-temp**2.))
                                endif
                            enddo
                        enddo
                    enddo
                    do nvar=1,23
                        do k=1,klev
                            do j=1,jm
                                do i=1,im
                                    var3d(i,j,k)=var3dav(i,j,k,nvar)
                                !   if(fsm2(i,j).eq.0.)var3d(i,j,k)=-1.e+10
                                enddo
                            enddo
                        enddo
                        write(*,*)'WRITE',nout,nvar
                    !    if(nvar.eq.1.or.nvar.eq.10.or.nvar.eq.17)then
                        call write_variable(ecoidfile2, var3d, aland)
                    !    if(nvar.eq.17)then
                    !        do j=1,jm
                    !            do i=1,im
                    !                var2d(i,j)=0.33*(var3d(i,j,1)+var3d(i,j,2)+var3d(i,j,3))
                    !            enddo
                    !        enddo
                    !        write(42)(var2d(i,1),i=1,im*jm)
                    !    endif
                        if(nvar == 20)then
                            do j=1,jm
                                do i=1,im
                                    var2d(i,j)=0.4*var3d(i,j,4)
                                enddo
                            enddo
                        !    write(42)(var2d(i,1),i=1,im*jm)
                            do j=1,jm
                                do i=1,im
                                    do k=1,klev
                                        if(var3d(i,j,k) < var3d(i,j,1)*0.01)then
                                            var2d(i,j)=0.5*(dpt(k-1)+dpt(k))
                                            goto 678
                                        endif
                                    enddo
                                    678 continue
                                    if(var2d(i,j) == 0.)var2d(i,j)=-1.E+10
                                enddo
                            enddo
                        !    write(42)(var2d(i,1),i=1,im*jm)
                        endif
                        1189 format(269(2X,e12.5))
                        do k=1,klev
                            do j=1,jm
                                do i=1,im
                                    var3dav(i,j,k,nvar)=0.
                                enddo
                            enddo
                        enddo
                    enddo
                    do nvar=1,1
                        do j=1,jm
                            do i=1,im
                                if(fsm(i,j) == 1.)then
                                    var2d(i,j)=var2dbotav(i,j,nvar)
                                else
                                    var2d(i,j)=aland
                                endif
                            enddo
                        enddo
                        call write_variable_hydro2d(ecoidfile2, var2d, aland)
                        do j=1,jm
                            do i=1,im
                                var2dbotav(i,j,nvar)=0.
                            enddo
                        enddo
                    enddo
                    call write_variable_hydro2d(ecoidfile2, var2dav, aland)
                    do j=1,jm
                        do i=1,im
                            var2dav(i,j)=0.
                        enddo
                    enddo
                    do nvar=1,5
                        do k=1,klev
                            do j=1,jm
                                do i=1,im
                                    var3d(i,j,k)=var3dhyav(i,j,k,nvar)
                                    if(nvar == 5)then
                                        var3d(i,j,k)=rho(i,j,k)
                                    endif
                            !        if(fsm2(i,j).eq.0.)var3d(i,j,k)=-1.e+10
                                enddo
                            enddo
                        enddo
                        if(nvar == 1 .OR. nvar == 2)then
                    !        write(43)(var3d(i,1,1),i=1,im*jm)
                        endif
                !        if(nvar.eq.3)then
                        call write_variable(ecoidfile2, var3d, aland)
                !        endif
                        do k=1,klev
                            do j=1,jm
                                do i=1,im
                                    var3dhyav(i,j,k,nvar)=0.
                                enddo
                            enddo
                        enddo
                    enddo
                !    write(43)(amldav(i,1),i=1,im*jm)
                    do j=1,jm
                        do i=1,im
                            amldav(i,j)=0.
                        enddo
                    enddo
                endif
            10 enddo
        5 enddo
    2 enddo
    close(33)
    close(44)
    stop
end program






subroutine read_variable(ifile, car, aland)
!! Read variable from file.
    include 'comblk.h'
    include 'fsm3d.h'

    real ::  Car(im,jm,klev)
    integer :: i, k, ifile

    do k=1, klev
        read(ifile) (Car(i,1,k),i=1,im*jm)
    enddo
    return
end subroutine read_variable







subroutine read_variable_hydro(ifile, car, aland)
!! Read hydrology variable from file.
    include 'comblk.h'
    include 'fsm3d.h'

    real ::  Car(im,jm,klev)
    integer :: i, k, ifile

    do k=1, klev
        read(ifile) (Car(i,1,k),i=1,im*jm)
    enddo
    return
end subroutine read_variable_hydro






subroutine read_variable_2d(ifile, sig, aland)
!! Read 2D variable from file.
    include 'comblk.h'
    include 'fsm3d.h'

    real ::  Sig(im,jm)
    integer :: i, k, ifile

    read(ifile) (Sig(i,1),i=1,im*jm)
    return
end subroutine read_variable_2d






subroutine write_variable(ifile, car, aland)
!! Write variable to file.
    include 'comblk.h'
    include 'fsm3d.h'

    real ::  Car(im,jm,klev)
    integer :: i, k, ifile

    do k=1, klev
        write(ifile) (Car(i,1,k),i=1,im*jm)
    enddo
    return
end subroutine write_variable






subroutine write_variable_hydro(ifile, car, aland)
!! Write hydrology variable to file.
    include 'comblk.h'
    include 'fsm3d.h'

    real ::   Car(im,jm,klev)
    integer :: i, k, ifile

    do k=1, klev
        write(ifile) (Car(i,1,k),i=1,im*jm)
    enddo
    return
end subroutine write_variable_hydro






subroutine write_variable_hydro2d(ifile, sig, aland)
!! Write 2D hydrology variable to file.
    include 'comblk.h'
    include 'fsm3d.h'

    real ::  Sig(im,jm)
    integer :: i, k, ifile

    write(ifile) (Sig(i,1),i=1,im*jm)
    return
end subroutine write_variable_hydro2d









subroutine DENS(SI,TI,RHOO,dpt)
!! Compute density.
!! T = potential temperature
    INCLUDE 'comblk.h'
    REAL*8 :: CR(IM,JM,KB),TR,SR,P
    REAL*8 :: RHOR(IM,JM,KB),RHO1,RHO2,RHO3,RHO4,RHO5
    dimension SI(IM,JM,KB),TI(IM,JM,KB),RHOO(IM,JM,KB)
    dimension DPT(KB)

!---compute density-1.025
    do 1 K=1,KBM1
        do 1 J=1,JM
            do 1 I=1,IM
        !        TR=TI(I,J,K)+10.E0
        !        SR=SI(I,J,K)+35.E0
        !---comvert decibars to bars
        !        P=-ZZ(K)*DT(I,J)*0.01*1.025*GRAV
                TR=TI(I,J,K)
                SR=SI(I,J,K)
                if(TR < 0.)goto 1
            !        P=DPT(K)*0.01*1.025*GRAV
                P=0.
                RHO1 = 999.842594 + 6.793952E-2*TR &
                - 9.095290E-3*TR**2 + 1.001685E-4*TR**3 &
                - 1.120083E-6*TR**4 + 6.536332E-9*TR**5
                RHO2 =   (0.824493 - 4.0899E-3*TR &
                + 7.6438E-5*TR**2 - 8.2467E-7*TR**3 &
                + 5.3875E-9*TR**4)
                RHO3 = RHO2*SR
                RHO4=  -5.72466E-3 + 1.0227E-4*TR - 1.6546E-6*TR**2
                RHO5 = RHO4*ABS(SR)**1.5  + 4.8314E-4 * SR**2
                RHOR(I,J,K)= RHO1 + RHO3 + RHO5
                CR(I,J,K) =1449.1+.0821*P+4.55*TR-.045*TR**2 &
                +1.34*(SR-35.)
                RHOR(I,J,K)=RHOR(I,J,K) + 1.E5*P/CR(I,J,K)**2 &
                *(1.-2.0*P/CR(I,J,K)**2)
                RHOO(I,J,K)=(RHOR(I,J,K)-1000.)
    1 enddo
    do 3 J=1,JM
        do 3 I=1,IM
            RHO(I,J,KB)=0.E0
    3 enddo
    return
end subroutine DENS







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

