program MAIN
!! This program makes files to be read by GRADS visualisation program. ??
    include 'sesame.fcm'
    include 'comblk.h'
    include 'ecooutputCO2.h'

    integer :: ecoidfile
    character ecofile*80
    dimension m_days_o(12),m_days_lp(12)
    real :: sst(im,jm)
    data m_days_o/31,28,31,30,31,30,31,31,30,31,30,31/
    data m_days_lp/31,29,31,30,31,30,31,31,30,31,30,31/
    character directory*100, dirday*11, outname*16, &
    outdirname*120,chara*1,outname2*16,indf1*2,indf2*1
    character(3) :: month(12)
    data month/'jan','feb','mar','apr','may','jun', &
    'jul','aug','sep','oct','nov','dec'/
    real :: aland
    data aland / -1.0E+10 /

    read(*,110) iday2,imo2,iyr22,directory
    write(*,*)iday2,imo2,iyr22,directory

!       directory='/data2/ktsiaras/POSEIDON_MED10_ECO_TESTRIV2/ &
!     MED_ERAdown/Run97hyd/'

!    10 format(1X,i2,1X,i2,1X,i4,1X,a80)
    110 format(1X,i2,1X,i2,1X,i4,1X,a80)
          
    ldf=0
    do ii=1,80
        chara=directory(ii:ii)
        if(chara /= ' ')then
            ldf=ldf+1
        endif
    enddo
    ihour2=18
    iwrite_at_end=0
    iyr2= iyr22 - 1900
    if (iyr22 >= 2000) iyr2 = iyr22 - 2000
    year=float(iyr22)
    if(mod(year,4.) == 0.)then
        mdays=m_days_lp(imo2)
    else
        mdays=m_days_o(imo2)
    endif
    i10day=1

!--read MED10 model grid
    open(550,FILE='../data/INITnew.DAT',FORM='UNFORMATTED')
    read(550) Z,ZZ,DZ,DZZ,ALON,ALAT,DX,DY,H
    close(550)
    do j=1,jm
        do i=1,im
            fsm(i,j)=1.
            if(h(i,j) <= 1.)fsm(i,j)=0.
        enddo
    enddo

!---read Ecology fields
    ilast=0
    icount=0
    i10=0
    call ECOZEROMM
    do 5 nouts=1,4000
!        if(iday2.le.mdays.and.iday2.gt.20)then
!            iaver=mdays-20
!        else
!            iaver=10
!        endif
!        if(iday2.eq.1)ilast=1
        iaver=mdays
!        iaver=8
!        if(imo2.eq.1.and.iday2.eq.4)ilast=1
        if(iday2 == 1)ilast=1
        imo2_p=imo2
        iyr2_p=iyr2
        call writeoutname(year,iyr2,imo2,iday2,ihour2, &
        outname,dirday,iwrite_at_end)

    !--open hydrodynamics
        outname(1:4)='OUTP'
        outdirname=directory(1:ldf)//dirday//outname
        write(*,*)outdirname
        open(3,file=outdirname,form='unformatted', &
        status='old',ERR=3000)
        do n=1,15
            read(3)
        enddo
        read(3)
        read(3)
        read(3)
        read(3) TIME,ELB,UA,VA,U,V,T,S,KH
        close(3)
        write(*,*)'TEST-ELB0',elemm(50,50)

    ! extract variables to plot from state vector
    ! 2(DO),3(CO2),5(phosphate),6(nitrate),7(ammonia),8(silicate),9(POC),13(DOC),16(DIAT),
    ! 20(nanoph),23(picopl),26(dynofl),29(mesozoo),30(microzoo),33(heter.dino),36(bact)
        icount=icount+1
        call HYDROACCUMMM
        write(*,*)'TEST-ELB1',elemm(50,50)

        goto 3100
        3000 continue
        if(icount == 0)then
            stop
        else
            ilast=1
        endif
        3100 continue

        if(icount == iaver .OR. ilast == 1)then
            aver=float(icount)

        !--open output file to write
            ecoidfile=33
            i10=i10+1
            write(indf1,'(i2.2)')iyr2_p
            write(indf2,'(i1)')4-i10
            ecofile=directory(1:ldf)//'hydro'//indf1//'.bin'
            write(*,*)ecofile
            open(33,file=ecofile,form='unformatted',access='append')
            ecofile=directory(1:ldf)//'ssh_sstMo'//indf1//'.bin'
            write(*,*)ecofile
            open(34,file=ecofile,form='unformatted',access='append')
            call HYDROFLUSHMM(aver,ecoidfile)
            write(*,*)'AVERAGING-ECO',icount,nouts,iaver
            call write_variable_hydro2d(34, elemm, aland)
            do j=1,jm
                do i=1,im
                    if(fsm(i,j) == 1.)then
                        sst(i,j)=tmm(i,j,1)
                    else
                        sst(i,j)=aland
                    endif
                enddo
            enddo
            call write_variable_hydro2d(34, sst, aland)
            call ECOZEROMM
            icount=0
        endif
        if(ilast == 1)stop
    5 enddo
end program

     








subroutine map1to3(im,jm,kb,fsm,N_COMP,ar1d, ar3d)
!! Map 1D variables from ERSEM to 3D ??
    real :: ar3d(im,jm,kb), ar1d(N_COMP)
    dimension fsm(im,jm)
    integer :: i,j,k,n

    iw = 0
    do k = 1, kb-1
        do j = 1, jm
            do i = 1, im
                if( fsm(i,j) > 0. ) then
                    iw = iw + 1
                    ar3d(i, j, k) = ar1d(iw)
                endif
            enddo
        enddo
    enddo

    if(iw /= N_COMP)then
        write(*,*)'ERROR in N_COMP',iw,N_COMP
        stop
    endif

    return
end subroutine map1to3









subroutine write_variable(ifile, var, aland)
!! Write variable to a file ??
!    implicit none
    include 'sesame_params.h'
    include 'comblk.h'
    include 'fsm3d.h'

    real :: var(N_COMP), Sig(im,jm,kb), Car(im,jm,klev), aland
    integer :: i, k, ifile

!   Zero the array
    do k = 1, kb
        do j = 1, jm
            do i = 1, im
                Sig(i, j, k) = val
            enddo
        enddo
    enddo

    call map1to3(im,jm,kb,fsm,N_COMP,var, Sig)
    call SIGMA2CART(Car, Sig, aland)
    do k=1, klev
        write(ifile) (Car(i,1,k),i=1,im*jm)
    enddo

    return
end subroutine write_variable








subroutine write_variable_hydro(ifile, sig, aland)
!! Write hydrology variable ??
    include 'comblk.h'
    include 'fsm3d.h'

    real ::  Sig(im,jm,kb), Car(im,jm,klev), aland
    integer :: i, k, ifile

    call SIGMA2CART(Car, Sig, aland)
    do k=1, klev
        write(ifile) (Car(i,1,k),i=1,im*jm)
    enddo

    return
end subroutine write_variable_hydro






subroutine write_variable_hydro2d(ifile, sig, aland)
!! Write 2D hydrology variable ??
    include 'comblk.h'
    include 'fsm3d.h'

    real ::  Sig(im,jm),  aland
    integer :: i, k, ifile
!    write(*,*)'TEST-ELB',sig(50,50),fsm(50,50)
           
    do j=1,jm
        do i=1,im
            if(fsm(i,j) == 0.)Sig(i,j)=aland
        enddo
    enddo
    write(ifile) (Sig(i,1),i=1,im*jm)
    return
end subroutine write_variable_hydro2d







          
subroutine writeoutname(year,iyr1,imo1,iday1,ihour1,outname,dirday,iwrite_at_end)
!! Write the name of a file using the date ??
    character dirday*11, outname*16
    dimension m_days_o(12),m_days_lp(12)
    data m_days_o/31,28,31,30,31,30,31,31,30,31,30,31/
    data m_days_lp/31,29,31,30,31,30,31,31,30,31,30,31/

    ihour=ihour1
    iday=iday1
    imo=imo1
    iyr=iyr1
    iyr_d=iyr1
    imo_d=imo1
    iyr_d=iyr1

    if(iwrite_at_end == 0)then
        iday_d=iday1
        goto 123
    endif
    iday_d=iday1-1
    if(iday_d == 0)then
        imo_d=imo_d-1
        if(imo_d == 0)then
            imo_d=12
            iyr_d=iyr_d-1
            if(iyr_d == -1)iyr_d=99
        endif
        if(mod(year,4.) == 0.)then
            iday_d=m_days_lp(imo_d)
        else
            iday_d=m_days_o(imo_d)
        endif
    endif

    123 continue

    write(*,*)'TEST',iyr1,imo1,iday1,ihour1

    write(outname,'(4x,3i2,1x,i2)') iday,imo,iyr,ihour
    outname(1:4)='ECOR'
    outname(11:11)='.'
    if (ihour <= 9) outname(12:12)='0'
    if (iday <= 9) outname(5:5)='0'
    if (imo <= 9) outname(7:7)='0'
    if (iyr <= 9) outname(9:9)='0'
    outname(14:16) = 'UTC'

    write(dirday,'(3i2)') iday_d,imo_d,iyr_d
    if (iday_d <= 9) dirday(1:1)='0'
    if (imo_d <= 9) dirday(3:3)='0'
    if (iyr_d <= 9) dirday(5:5)='0'
    dirday(7:11) = '_OCE/'

    iday1=iday1-1
    if(iday1 == 0)then
        imo1=imo1-1
        if(imo1 == 0)then
            imo1=12
            iyr1=iyr1-1
            if(iyr1 == -1)iyr1=99
        endif
        if(mod(year,4.) == 0.)then
            iday1=m_days_lp(imo1)
        else
            iday1=m_days_o(imo1)
        endif
    endif
    return
end subroutine writeoutname







subroutine ECOZEROMM
!! Initialize all accumulators to zero.
!    implicit none
    include 'comblk.h'
    include 'sesame_params.h'
    include 'ecooutputCO2.h'

    integer :: n

    do 10 n = 1, N_COMP
        n1pmm(n) = 0.0
        n3nmm(n) = 0.0
        n4nmm(n) = 0.0
        n5smm(n) = 0.0
        p1cmm(n) = 0.0
        p2cmm(n) = 0.0
        p3cmm(n) = 0.0
        p4cmm(n) = 0.0
        r1cmm(n) = 0.0
        z4cmm(n) = 0.0
        z5cmm(n) = 0.0
        z6cmm(n) = 0.0
        b1cmm(n) = 0.0
        r6cmm(n) = 0.0
        o2omm(n) = 0.0
        chlmm(n) = 0.0
        netppmm(n) = 0.0
        prodB1mm(n) = 0.0
        xepsmm(n) = 0.0
        eirmm(n) = 0.0
    10 enddo

    do 20 i=1,im
        do 20 j=1,jm
            elemm(i,j) = 0.0
            do 30 k=1,kb
                umm(i,j,k) = 0.0
                vmm(i,j,k) = 0.0
                tmm(i,j,k) = 0.0
                smm(i,j,k) = 0.0
                akhmm(i,j,k) = 0.0
            30 enddo
    20 enddo
    return
end subroutine ECOZEROMM








subroutine ECOACCUMMM
!! Accumulate monthly means of ecology variables.
!    implicit none
    include 'comblk.h'
    include 'sesame.fcm'
    include 'ecooutputCO2.h'

    integer :: n

!---Accumulate monthly means:
    do 10 n = 1, N_COMP
        n1pmm(n) = n1pmm(n) + stvar(n1p,n)
        n3nmm(n) = n3nmm(n) + stvar(n3n,n)
        n4nmm(n) = n4nmm(n) + stvar(n4n,n)
        n5smm(n) = n5smm(n) + stvar(n5s,n)
        p1cmm(n) = p1cmm(n) + stvar(p1c,n)
        p2cmm(n) = p2cmm(n) + stvar(p2c,n)
        p3cmm(n) = p3cmm(n) + stvar(p3c,n)
        p4cmm(n) = p4cmm(n) + stvar(p4c,n)
        r1cmm(n) = r1cmm(n) + stvar(r1c,n)
        z4cmm(n) = z4cmm(n) + stvar(z4c,n)
        z5cmm(n) = z5cmm(n) + stvar(z5c,n)
        z6cmm(n) = z6cmm(n) + stvar(z6c,n)
        b1cmm(n) = b1cmm(n) + stvar(b1c,n)
        r6cmm(n) = r6cmm(n) + stvar(r6c,n)
        o2omm(n) = o2omm(n) + stvar(o2o,n)
        chlmm(n) = chlmm(n) + chl(n)
        netppmm(n) = netppmm(n) + netpp(n)
        prodB1mm(n) = prodB1mm(n) + prodB1(n)
        xepsmm(n) = xepsmm(n) + xeps(n)
        eirmm(n) = eirmm(n) + eir(n)
    10 enddo

    do 20 i=1,im
        do 20 j=1,jm
            if(fsm(i,j) /= 1.) goto 20
            elemm(i,j) = elemm(i,j) + elb(i,j)
            do 30 k=1,kb
                umm(i,j,k) = umm(i,j,k) + u(i,j,k)
                vmm(i,j,k) = vmm(i,j,k) + v(i,j,k)
                tmm(i,j,k) = tmm(i,j,k) + t(i,j,k)+10.
                smm(i,j,k) = smm(i,j,k) + s(i,j,k)+35.
                akhmm(i,j,k) = akhmm(i,j,k) + kh(i,j,k)
            30 enddo
    20 enddo
            
    write(*,*)'TEST-ELB2',elemm(50,50)
    return
end subroutine ECOACCUMMM







subroutine HYDROACCUMMM
!! Accumulate monthly means of hydrology variables.
!    implicit none
    include 'comblk.h'
    include 'sesame.fcm'
    include 'ecooutputCO2.h'

    integer :: n

    do 20 i=1,im
        do 20 j=1,jm
            if(fsm(i,j) /= 1.) goto 20
            elemm(i,j) = elemm(i,j) + elb(i,j)
            if(elb(i,j) > 0.5)write(*,*)'HIGH',i,j,h(i,j),elb(i,j)
            do 30 k=1,kb
                umm(i,j,k) = umm(i,j,k) + u(i,j,k)
                vmm(i,j,k) = vmm(i,j,k) + v(i,j,k)
                tmm(i,j,k) = tmm(i,j,k) + t(i,j,k)+10.
                smm(i,j,k) = smm(i,j,k) + s(i,j,k)+35.
                akhmm(i,j,k) = akhmm(i,j,k) + kh(i,j,k)
            30 enddo
    20 enddo

    write(*,*)'TEST-ELB2',elemm(50,50),elb(50,50)
    return
end subroutine HYDROACCUMMM








subroutine ECOFLUSHMM(aver,ecoidfile)
!! Normalize ecology variables and write them. ??
!      implicit none
    include 'comblk.h'
    include 'sesame_params.h'
    include 'ecooutputCO2.h'

    real :: aland
    data aland / -1.0E+10 /
    integer :: n

!---Normalize Ecology Variables
    do 10 n = 1, N_COMP
        n1pmm(n) = n1pmm(n)/aver
        n3nmm(n) = n3nmm(n)/aver
        n4nmm(n) = n4nmm(n)/aver
        n5smm(n) = n5smm(n)/aver
        p1cmm(n) = p1cmm(n)/aver
        p2cmm(n) = p2cmm(n)/aver
        p3cmm(n) = p3cmm(n)/aver
        p4cmm(n) = p4cmm(n)/aver
        r1cmm(n) = r1cmm(n)/aver
        z4cmm(n) = z4cmm(n)/aver
        z5cmm(n) = z5cmm(n)/aver
        z6cmm(n) = z6cmm(n)/aver
        b1cmm(n) = b1cmm(n)/aver
        r6cmm(n) = r6cmm(n)/aver
        o2omm(n) = o2omm(n)/aver
        chlmm(n) = chlmm(n)/aver
        netppmm(n) = netppmm(n)/aver
        prodB1mm(n) = prodB1mm(n)/aver
        xepsmm(n) = xepsmm(n)/aver
        eirmm(n) = eirmm(n)/aver
    10 enddo
    write(*,*)'TEST-ELB',elemm(50,50)

    do 20 i=1,im
        do 20 j=1,jm
            if(fsm(i,j) /= 1.) goto 20
            elemm(i,j) = elemm(i,j)/aver
            do 30 k=1,kb
                umm(i,j,k) = umm(i,j,k)/aver
                vmm(i,j,k) = vmm(i,j,k)/aver
                tmm(i,j,k) = tmm(i,j,k)/aver
                smm(i,j,k) = smm(i,j,k)/aver
                akhmm(i,j,k) = akhmm(i,j,k)/aver
            30 enddo
    20 enddo

    write(*,*)'TEST-ELBB',elemm(50,50)

!---Write Ecology Variables
    call write_variable(ecoidfile, n1pmm, aland)
    call write_variable(ecoidfile, n3nmm, aland)
    call write_variable(ecoidfile, n4nmm, aland)
    call write_variable(ecoidfile, n5smm, aland)
    call write_variable(ecoidfile, p1cmm, aland)
    call write_variable(ecoidfile, p2cmm, aland)
    call write_variable(ecoidfile, p3cmm, aland)
    call write_variable(ecoidfile, p4cmm, aland)
    call write_variable(ecoidfile, r1cmm, aland)
    call write_variable(ecoidfile, z4cmm, aland)
    call write_variable(ecoidfile, z5cmm, aland)
    call write_variable(ecoidfile, z6cmm, aland)
    call write_variable(ecoidfile, b1cmm, aland)
    call write_variable(ecoidfile, r6cmm, aland)
    call write_variable(ecoidfile, o2omm, aland)
    call write_variable(ecoidfile, chlmm, aland)
    call write_variable(ecoidfile, netppmm, aland)
    call write_variable(ecoidfile, prodB1mm, aland)
    call write_variable(ecoidfile, xepsmm, aland)
    call write_variable(ecoidfile, eirmm, aland)

    call write_variable_hydro2d(ecoidfile, elemm, aland)
    call write_variable_hydro(ecoidfile, umm, aland)
    call write_variable_hydro(ecoidfile, vmm, aland)
    call write_variable_hydro(ecoidfile, tmm, aland)
    call write_variable_hydro(ecoidfile, smm, aland)
    call write_variable_hydro(ecoidfile, akhmm, aland)

    return
end subroutine ECOFLUSHMM








subroutine HYDROFLUSHMM(aver,ecoidfile)
!! Normalize hydrology variables and write them. ??
!    implicit none
    include 'comblk.h'
    include 'sesame_params.h'
    include 'ecooutputCO2.h'

    real :: aland
    data aland / -1.0E+10 /
    integer :: n

    do 20 i=1,im
        do 20 j=1,jm
            if(fsm(i,j) /= 1.) goto 20
            elemm(i,j) = elemm(i,j)/aver
            do 30 k=1,kb
                umm(i,j,k) = umm(i,j,k)/aver
                vmm(i,j,k) = vmm(i,j,k)/aver
                tmm(i,j,k) = tmm(i,j,k)/aver
                smm(i,j,k) = smm(i,j,k)/aver
                akhmm(i,j,k) = akhmm(i,j,k)/aver
            30 enddo
    20 enddo
    call write_variable_hydro2d(ecoidfile, elemm, aland)
    call write_variable_hydro(ecoidfile, umm, aland)
    call write_variable_hydro(ecoidfile, vmm, aland)
    call write_variable_hydro(ecoidfile, tmm, aland)
    call write_variable_hydro(ecoidfile, smm, aland)
    call write_variable_hydro(ecoidfile, akhmm, aland)

    return
end subroutine HYDROFLUSHMM







subroutine TRANSLATE(TITLE,IYR_N,IMO_N,IDAY_N)
!! ??
    character(10) :: TITLE,SYM(0:9)*1
    character(1) :: AUX
    dimension IHOLD(3)
    integer :: QUANT(4)
    data QUANT/1,10,100,1000/
    data SYM/'0','1','2','3','4','5','6','7','8','9'/

    K = 1
    L = 1

    do N = 1,3
        IHOLD(N) = 0
    enddo
    do 200 I = 10,1,-1
        AUX = TITLE(I:I)
        if (AUX == '/') then
            L = 1
            K = K+1
            goto 200
        endif
        do M = 0,9
            if (AUX == SYM(M)) then
                MULT = M
                goto 100
            endif
        enddo
        100 continue
        IHOLD(K) = IHOLD(K) + MULT*QUANT(L)
        L = L + 1
    200 enddo

    IYR_N = IHOLD(1)
    IMO_N = IHOLD(2)
    IDAY_N = IHOLD(3)
    return
end subroutine TRANSLATE

