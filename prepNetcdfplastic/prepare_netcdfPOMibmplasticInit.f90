program prepare_netcdf_ibmPlastic
!!  This program reads prepares the initial ibm file for plastics 
!!  number of particles for each size class is based on a homogeanous background concentration 
!!  equal to the median of available in situ data (see plasticparams.dat)
    use poseidon_common_subroutines
    use ibm_common_subroutines
    include "netcdf.inc"
    include 'grid.h'
    include 'sesame_params.h'
    include 'ibmall.fcm'
    parameter  (ip=0)                               !add more particles in the same grid box (e.g ip=1->np=9)
    parameter  (nadd=2*ip+1)           
    parameter  (iskip=1)                            !particles in every iskip grid boxes
    parameter (unitArea10=9.6450272E+01)            !grid box area (km2)
    real :: xilon(n_upper$g),yilat(n_upper$g)
    integer :: ch
    integer :: nsi_id,species_id,nfi_id,variant_id,xsi_id,ysi_id, &
    ref_id,iarea_id,fdep_id,icountdays_id,rbt_id,fou_id,dista_id

    integer :: alon_id,alat_id,h_id,fsm_id

    real :: countvariant(6,2)
    real :: concbiom(im,jm),conc(im,jm,nvariantsmax,nspecies),conc2d(im,jm)
    parameter (MINCHK=-900)
    integer :: NCID1
    character CFILE1*128

!---read model grid (used in subroutine init_lagrSI)----------------
!---h=bathymetry, alon/alat=longitude,latitude, fsm=sea mask(=0 for land)
   status = NF_OPEN('../in/pom.gridplastic.nc',0, ncid)
    write(*,*)'OK',filename,ncid
!---get variables
    status=nf_inq_varid(ncid,'h',h_id)
    status=nf_inq_varid(ncid,'longitude',alon_id)
    status=nf_inq_varid(ncid,'latiitude',alat_id)
    status=nf_inq_varid(ncid,'fsm',fsm_id) 
    status=nf_get_var_real(ncid,h_id,h)
    status=nf_get_var_real(ncid,alon_id,alon)
    status=nf_get_var_real(ncid,alat_id,alat)
    status=nf_get_var_real(ncid,fsm_id,fsm)
    status = NF_CLOSE (NCID)

    


!---set uniform SI's locations
    xstep=xres
    ystep=yres
    if(ip > 0)then
        xadd=xstep/(4*ip)
        yadd=ystep/(4*ip)
    else
        xadd=0.
        yadd=0.
    endif
    if(iskip > 1)then
        xskip=xstep*0.5*(iskip-1)
        yskip=ystep*0.5*(iskip-1)
    else
        xskip=0.
        yskip=0.
    endif
    xadd=0.
    yadd=0.
    xskip=0.
    yskip=0.
    iw = 0
    do j = 1, jm
        do i = 1, im
            if( fsm(i,j) > 0.)then
                if(mod(i-1,iskip) == 0. .AND. mod(j-1,iskip) == 0)then
                    do ii=-ip,ip
                        do jj=-ip,ip
                            iw = iw + 1
                            xilon(iw) = alon(i,j) +ii*xadd + xskip
                            yilat(iw) = alat(i,j) +jj*yadd + yskip
                    !        write(*,*)i,j,iw,xilon(iw),yilat(iw),xilon(1),yilat(1)
                        enddo
                    enddo
                endif
            endif
        enddo
    enddo
    nsea=iw

!---find coast
    do j=1,jm
        do i=1,im
            coast(i,j)=0
            if(fsm(i,j) == 0.)then
                nav=0
                do jj=-1,1
                    do ii=-1,1
                        if(i+ii > 0 .AND. i+ii <= im .AND. j+jj > 0 &
                         .AND. j+jj <= jm)then
                            nav=(fsm(i+ii,j+jj))+nav
                        endif
                    enddo
                enddo
                if (nav > 0.)then
                    coast(i,j)=1
                endif
            endif
        enddo
    enddo

    ncoast=iw-nsea

!    nareas=ncoast+nsea

    nareas=nsea
    write(*,*)'NAREAS',nareas,nsea,ncoast

    call initfish(nareas,init_nsi)
    write(*,*)'initfish',nareas,xilon(1),yilat(1)

    call init_lagr(xilon,yilat,nareas,nsea,init_nsi)
    write(*,*)'initlagr',nsi,nareas,xilon(1),yilat(1)
        
    do n=1,nsi
        ivariant=variant(n)
        nsp=species(n)
        if(ivariant > 0)then
            countvariant(ivariant,nsp)=countvariant(ivariant,nsp)+1
        endif
    enddo
    do nsp=1,2
        do ivariant=1,6
            write(*,*)'VARIANTS',nsp,ivariant,countvariant(ivariant,nsp)
        enddo
    enddo
    xs=alon(2,1)-alon(1,1)
    ic=0
    do ch=1,nsi
        ivariant=variant(ch)
        nsp=species(ch)
        ii=nint((xsi(ch)-xleft-0.5*xs)/xs)+1
        jj=nint((ysi(ch)-ybot-0.5*xs)/xs)+1
        if(variant(ch) == 1 .AND. species(ch) == 1)then
            concbiom(ii,jj)=concbiom(ii,jj)+ &
            nfi(ch)*xfi(ch)*1.e-6/unitArea10
        endif
        if(species(ch) /= 0)then
            conc(ii,jj,ivariant,nsp)=conc(ii,jj,ivariant,nsp)+ &
            nfi(ch)/unitArea10
        endif
    enddo

    call PRXY('BIOM',TIME,concbiom,im,1,jm,1,1.E-4)

    open(12,file='testibm.bin',form='unformatted')
    do nsp=1,nspecies
        do ivariant=1,nvariants(nsp)
            do j=1,jm
                do i=1,im
                    if(conc(i,j,ivariant,nsp) == 0.)conc(i,j,ivariant,nsp)=-1.E+10
                    conc2d(i,j)=conc(i,j,ivariant,nsp)
                enddo
            enddo
            if(nsp == 2 .AND. ivariant == 4)then
                call PRXY('BIOM',TIME,conc2d,im,1,jm,1,1.E-0)
            endif
            write(12)(conc(i,1,ivariant,nsp),i=1,im*jm)
        enddo
    enddo
    close(12)

    nsi=nmax
    write(*,*)'NSI=',nsi
!    do n=1,1000
!        write(*,*)n,species(n),xsi(n),xfi(n),variant(n)
!    enddo
!---------------
    icount=0
    do k=1,nsi
        if(ref(k) == 99)then
            icount=icount+1
            isi(k)=int((xsi(k)-xleft)/xres)+1
            jsi(k)=int((ysi(k)-ybot)/xres)+1
            if(fsm(isi(k),jsi(k)) == 1)write(*,*)'TEST-COAST', &
            k,xsi(k),ysi(k)
        endif
    enddo
    write(*,*)'REF99',icount

!    do k=1,nsi
!        variant(k)=0
!    enddo

!---open NetCDF files
!    CFILE1='../data/pom.ibmplastic0.nc'
    CFILE1='../data/pom.ibmplastic.nc'
    istatus=NF_CREATE(CFILE1(1:(INDEX(CFILE1,' '))-1),NF_SHARE,NCID1)!create dataset
    istatus = NF_DEF_DIM(NCID1, 'nsi',nsi, nsi_id)
    istatus = NF_REDEF(NCID1)
    istatus = NF_DEF_VAR(NCID1,'species',NF_INT,1,nsi_id,species_id)
    istatus = NF_DEF_VAR(NCID1,'nfi',NF_FLOAT,1,nsi_id,nfi_id)
    istatus = NF_DEF_VAR(NCID1,'variant',NF_INT,1,nsi_id,variant_id)
    istatus = NF_DEF_VAR(NCID1,'xsi',NF_FLOAT,1,nsi_id,xsi_id)
    istatus = NF_DEF_VAR(NCID1,'ysi',NF_FLOAT,1,nsi_id,ysi_id)
    istatus = NF_DEF_VAR(NCID1,'fdep',NF_FLOAT,1,nsi_id,fdep_id)
    istatus = NF_DEF_VAR(NCID1,'ref',NF_INT,1,nsi_id,ref_id)
    istatus = NF_DEF_VAR(NCID1,'iarea',NF_INT,1,nsi_id,iarea_id)
    istatus = NF_DEF_VAR(NCID1,'icountdays',NF_INT,1,nsi_id, &
    icountdays_id)
    istatus = NF_DEF_VAR(NCID1,'dista',NF_FLOAT,1,nsi_id, &
    dista_id)
    istatus = NF_DEF_VAR(NCID1,'rbt',NF_FLOAT,1,nsi_id, &
    rbt_id)
    istatus = NF_DEF_VAR(NCID1,'fou',NF_FLOAT,1,nsi_id, &
    fou_id)
    istatus = NF_ENDDEF(ncid1, iret)
    istatus = nf_put_var_int(ncid1,species_id,species)
    istatus = nf_put_var_real(ncid1,nfi_id,nfi)
    istatus = nf_put_var_int(ncid1,variant_id,variant)
    istatus = nf_put_var_real(ncid1,xsi_id,xsi)
    istatus = nf_put_var_real(ncid1,ysi_id,ysi)
    istatus = nf_put_var_real(ncid1,fdep_id,fdep)
    istatus = nf_put_var_int(ncid1,ref_id,ref)
    istatus = nf_put_var_int(ncid1,iarea_id,iarea)
    istatus = nf_put_var_int(ncid1,icountdays_id,icountdays)
    istatus = nf_put_var_real(ncid1,rbt_id,rbt)
    istatus = nf_put_var_real(ncid1,fou_id,fou)
    istatus = nf_put_var_real(ncid1,dista_id,dista)
    istatus = NF_CLOSE (NCID1)

    icount=0
    do n=1,nsi
        if(variant(n) == 1 .AND. species(n) == 1)then
            icount=icount+1
            write(*,*)'TESTINI',n,icount,ysi(n)
        endif
    enddo
    stop
end program

subroutine initfish(nareas,init_nsi)
!! Set initial Super Individual populations,
    use poseidon_common_subroutines
    use ibm_common_subroutines
    include 'ibmall.fcm'
    include 'init.h'
    integer :: ch,i,nar,nn
    integer :: icountempty
    integer :: ipat,jpat,ii,init_nsi(nspecies),ic, &
    init_nsiall
    integer :: nsp
!---set initial SI populations

    init_nsiall=0
    write(*,*)'Initialise plastic, SPECIES='
    do nsp=1,nspecies
        init_nsi(nsp)=0
        write(*,*)'SPECIES=',nspecies,nvariants(1),nvariants(2)
        ic=0
        do i=1,nvariants(nsp)
            write(*,*)'VARIANT',i,nsiInit(i,nsp)*nareas
            init_nsi(nsp)=init_nsi(nsp)+nsiInit(i,nsp)*nareas
            init_nsiall=init_nsiall+nsiInit(i,nsp)*nareas
            if(nsiInit(i,nsp) > 0)then
                ic=ic+1
                icvariant(i,nsp)=ic
            endif
        enddo
        write(*,*)'NSIs',init_nsi(nsp)
    enddo
    write(*,*)'NSIs',init_nsiall

!    nsi=init_nsiall
    nsi=nmax
    do n=1,nsi
        variant(n)=0
    enddo
    icount=0
    do n=1,nsi
        if(variant(n) == 2)icount=icount+1
    enddo
    write(*,*)'TEST1',icount,nmax

!---Initialize ibm's of different classes
    iadd=0
    icount=0
    do nsp=1,nspecies
        write(*,*)'Species',nsp
        do i=1,nvariants(nsp)
        !   do i=1,1
            nvariantloc=nvariants(nsp)
            write(*,*)'Variant',i,icvariant(i,nsp),nsiInit(i,nsp)
            do nn=1,nsiInit(i,nsp)*nareas
                ch=(nvariantloc*(nn-1)+icvariant(i,nsp)+iadd)
                variant(ch)=i
                species(ch)=nsp
                xfi(ch)=xfiInit(i,nsp)
                nfi(ch)=nfiInit(i,nsp)
                icountdays(ch)=0
                dista(ch)=0.
                if(variant(ch) <= 3 .AND. species(ch) == 2)then
                    fou(ch)=dens(i,nsp)
                else
                    fou(ch)=0.
                endif
!                write(*,*)'Variant',i,ch,real(xfi(ch)),real(nfi(ch))
                write(*,*)'Variant',i,ch,xfi(ch),nfi(ch)
                if(variant(ch) == 2)icount=icount+1
            enddo
        enddo
        iadd=iadd+init_nsi(nsp)
    enddo
    write(*,*)'TEST2',icount,nmax

    icount=0
    do n=1,nsi
        if(variant(n) == 2)icount=icount+1
    enddo
    write(*,*)'TEST3',icount,nmax

    do n=1,nmax
        nempty(n)=0
    enddo
    icountempty=0
    do n=1,nsi
        if(variant(n) == 0)then
            icountempty=icountempty+1
        endif
    enddo
    write(*,*)icountempty,nmax

    call totfish(1)
    return
end subroutine initfish

    


subroutine totfish(iwrite)
!! write total number of SIs(tsivariant), particles(tnvariant) and weight(tbvariant) for each class
!! calculate empty SIs 
    use poseidon_common_subroutines
    use ibm_common_subroutines
    include 'ibmall.fcm'
    integer :: ch,i,nn,iwres
    integer :: icountempty,iwrite

    do nsp=1,nspecies
        do i=1,nvariants(nsp)
            weivariant(i,nsp)=0.
            tbvariant(i,nsp)=0.
            tnvariant(i,nsp)=0.
            tsivariant(i,nsp)=0.
        enddo
    enddo

    do n=1,nmax
        nempty(n)=0
    enddo
    icountempty=0
    do n=1,nsi
        if(variant(n) == 0)then
            icountempty=icountempty+1
            nempty(icountempty)=n
        else
            nsp=species(n)
            ivariant=variant(n)
            tbvariant(ivariant,nsp)=tbvariant(ivariant,nsp)+nfi(n)*xfi(n)
            tnvariant(ivariant,nsp)=tnvariant(ivariant,nsp)+nfi(n)
            tsivariant(ivariant,nsp)=tsivariant(ivariant,nsp)+1
            weivariant(ivariant,nsp)=weivariant(ivariant,nsp)+nfi(n)*xfi(n)
        endif
    enddo

    do n=1,nspecies
        do i=1,nvariants(n)
            if(tnvariant(i,n) > 0.)then
                weivariant(i,n)=weivariant(i,n)/tnvariant(i,n)
            else
                weivariant(i,n)=0.
            endif
        enddo
    enddo

    if(iwrite == 1)then
        write(*,*)'NSI=',nsi
        write(*,*)'EMPTY=',icountempty
        do n=1,nspecies
            write(*,*)'SPECIES',n
            adultbiom=0.
            do ivariant=1,nvariants(n)
                write(*,*)'TOTAL BIOMASS(tn)',ivariant,tbvariant(ivariant,n)*1.e-6
                write(*,*)'PARTICLES',ivariant,tnvariant(ivariant,n)
                write(*,*)'SIs',ivariant,tsivariant(ivariant,n)
                write(*,*)'PARTICLE WEIGHT(gr)',ivariant,weivariant(ivariant,n)
            enddo
        enddo
    endif
    return
end subroutine totfish


subroutine init_lagr(xilon,yilat,nareas,nsea,init_nsi)
!! Initialize lagrangian particles
    use poseidon_common_subroutines
    use ibm_common_subroutines
    include 'ibmall.fcm'
    include 'grid.h'
    include 'lagr.h'
    include 'init.h'
    real :: xilon(nareas),yilat(nareas)
    real :: xi(nareas),yi(nareas) !position in cartesian coordinates (m)
    real :: xcg(nareas),ycg(nareas),vcg(nareas) !centre of gravity
    real :: xxsi,yysi
    integer :: iwres,n1,n2,mb,nvariantloc
    integer :: init_nsi(nspecies)
    real ::        pi,                g
    parameter ( pi = 3.1415926536, g = 9.81 )

    do k=1,nsi
        ref(k)=-9
    end do

    do ms=1,nareas

!---Proceed with Initialization - Get Geographic coordinates
!            and Set Cartesian and Grid coordinates

!---xilon(ms),yilat(ms) specified in ibm.f90
!        print *, 'Source=',ms
!        print *, 'Initial Geographic coordinates=', xilon(ms),yilat(ms)

!---Set CARTESIAN coordinates (meters)
        call Geo2Cart(xi(ms), yi(ms), alon(1,1)-(stepx/2.), &
            alat(1,1)-(stepy/2.), xilon(ms),yilat(ms), pi, dg)
!        write(*,*)ms,xilon(ms),yilat(ms),xi(ms),yi(ms)
    enddo   !end do ms=1,nareas

!---Proceed with intializations
! Before other initializations find different bc on land boxes
    call land_bc(lbc, h, im, jm, 1.)
!---assign retention time for coastal areas
    do j = 1, jm
        do i = 1, im
            rbc(i, j) = 0.
            if(h(i, j) <= 1.) rbc(i, j) = rtime(1,1)
        end do
    end do
    write(*,*)'OK2',nareas,xilon(1),yilat(1)
!---position of each particle
    iadd=0
    k=0
    do 5 nsp=1,nspecies
        write(*,*)'SPECIES',nsp
        do 10 ivariant=1,nvariants(nsp)
            write(*,*)'VARIANT',ivariant
            nvariantloc=nvariants(nsp)
            do 20 ms=1,nsea
    !            write(*,*)'AREA',ms
                k=(nvariantloc*(ms-1)+icvariant(ivariant,nsp)+iadd)
                write(*,*)'AREA',ms,nvariantloc,k
                if(variant(k) == 0)goto 122
    !            iarea(k)=ms
                iarea(k)=1
                50 continue
                r1 = rand(0)
                r2 = rand(0)
                xn(k) = xi(ms)
                yn(k) = yi(ms)
!                xn(k) = xi(ms) + r1 * rin * cos(2. * pi * r2)
!                yn(k) = yi(ms) + r1 * rin * sin(2. * pi * r2)
                zn(k) = 10.
!                xn(k) = xi(ms)
!                yn(k) = yi(ms)
!                zn(k) = 10.
                xo(k) = xn(k)
                yo(k) = yn(k)
                zo(k) = zn(k)
                fdep(k)=0.
                call FCoord(io, jo, alon(1,1)-(stepx/2.), &
                    alat(1,1)-(stepy/2.), stepx, stepy, &
                    xo(k), yo(k), pi, dg)
                in = io
                jn = jo
    !            write(*,*)'IBM',k,xilon(ms),yilat(ms),io,jo,h(io, jo)
                if(h(io, jo) <= 1.) goto 50
                if(h(io, jo) > 1.) then
                    ref(k) = 0
                    rbt(k) = 0.
                end if
            !---initialize ibm (i,j) and position
                call Cart2Geo(xsi(k), ysi(k), alon(1,1)-(stepx/2.), &
                alat(1,1)-(stepy/2.), xn(k), yn(k), pi, dg)
            !    write(*,*)'XSI1',k,xsi(k),ysi(k),ysi(100390)
                call Geo2Cart(xn(k), yn(k), alon(1,1)-(stepx/2.), &
                alat(1,1)-(stepy/2.), xsi(k),ysi(k), pi, dg)
                write(*,*)'XSI',k,xsi(k),ysi(k)
            !---Initialization of the gravity centers
                xcg(ms)=0.
                ycg(ms)=0.
                vcg(ms)=0.
                122 continue
                continue     !ms=1,nareas / Position of each particle
            20 enddo
        10 enddo
        iadd=iadd+init_nsi(nsp)
    5 enddo
!    write(*,*)'TESTXSI2',ysi(100390)

!    do k=1,nsi
!        call Geo2Cart(xo(k), yo(k), alon(1,1)-(stepx/2.),alat(1,1)-(stepy/2.),xsi(k),ysi(k),pi,dg)
!        write(*,*)'IBM RES',species(k),xsi(k),ysi(k)
!        zo(k)=fdep(k)
!        if(variant(k).gt.0.and.ref(k).ne.0)write(*,*)'INIT-LAGR',k,xsi(k),ysi(k),variant(k),ref(k)
!    enddo
!    write(*,*)'TESTXSI',k,nsp,ivariant,ysi(100390)
    return
    123 continue
!    write(*,*)'INITIALIZE IBM from RESTART'
!    do k=1,nsi
!        call Geo2Cart(xo(k), yo(k), alon(1,1)-(stepx/2.),alat(1,1)-(stepy/2.),xsi(k),ysi(k),pi,dg)
!        write(*,*)'IBM RES',xsi(k),ysi(k)
!        zo(k)=fdep(k)
!    enddo
    return
end subroutine init_lagr


