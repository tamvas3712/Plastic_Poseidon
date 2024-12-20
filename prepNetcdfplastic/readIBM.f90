program read_IBM
    include 'ibmall.fcm'
!    include 'grid.h'
    parameter (im=431,jm=156,kb=25)
    parameter (nobsmicro=683,nobsmacro=475)
    parameter (xleft=-7,ybot=30.25,xs=0.1,ys=0.1) !model grid
    include 'fsm3d.h'
    include 'lagr.h'
    parameter (im2=im,jm2=jm)                            !downscaled grid for better resolution
    parameter (unitArea10=9.6450272E+01)                   !km2
    parameter (unitArea20=unitArea10/4.)
    parameter (unitArea240=unitArea10/576.)
    parameter (unitArea=unitArea10)
    integer :: isiind(nmax)
    integer :: nsi_id,species_id,xfi_id,nfi_id,variant_id,xsi_id,ysi_id, &
    iflagVariant_id,weirate_id,rbufcum_id, &
    conxfi_id,ref_id,iarea_id,indsi_id,rbt_id,fou_id,nschoo_id, &
    icountEgg_id,keepeggs_id,adday_coho_id,fdep_id,icountdays_id
!---fish SI read
    integer :: ch
    real :: h(im,jm),alon(im,jm),alat(im,jm),dx(im,jm),dy(im,jm), &
    z(kb),zz(kb),dz(kb),dzz(kb)
!---2-D distributions
    real :: concbiom(im,jm,nvariantsmax,nspecies),     & !biomass gr/Km2
    conc(im,jm,nvariantsmax,nspecies),         & !particles/Km2
    conc1(im,jm,nvariantsmax,nspecies),         & !particles/Km2
    conc2(im,jm,nvariantsmax,nspecies),         & !particles/Km2
    concbeach(im,jm,nvariantsmax,nspecies),         & !particles/Km2
    concbottom(im,jm,nvariantsmax,nspecies),         & !particles/Km2
    sidep(im,jm,nvariantsmax,nspecies),concbt(im,jm), &
    countsi1(im,jm,nvariantsmax,nspecies), &
    countsi2(im,jm,nvariantsmax,nspecies)
    real :: sink2d(im,jm,nvariantsmax,nspecies), &
    density(im,jm,nvariantsmax,nspecies),fou2d(im,jm,nvariantsmax,nspecies), &
    countsink(im,jm,nvariantsmax,nspecies), &
    sidays(im,jm,nvariantsmax,nspecies),sifou(im,jm,nvariantsmax,nspecies), &
    sidist(im,jm,nvariantsmax,nspecies)
    real :: conc3dav(37,nvariantsmax),conc3d(im,jm,37,nvariantsmax)
    character directory*140,fname1*26,fname*160,fname2*170,chara*1
    character fishname*80
    real :: coast(im,jm),fsm(im,jm)
    real :: fdepav1(365),fdepav2(365),wbav1(365),wbav2(365), &
    hfouav1(365),hfouav2(365),denav1(365),denav2(365), &
    count1(365),count2(365)
    real :: concmean(nvariantsmax,nspecies),concmeanBea(nvariantsmax,nspecies)
    real :: totbiom(nspecies),totbiom100(nspecies),totbiomBea(nspecies), &
    totbiomBot(nspecies)
!---model values at data points
    real :: conc_model_macro(5),conc_model_micro(2)

    open(50,FILE='../data/INITnew.DAT',FORM='UNFORMATTED')
    read(50) Z,ZZ,DZ,DZZ,ALON,ALAT,DX,DY,H
    close(50)

    do j=1,jm
        do i=1,im
            fsm(i,j)=1.
            if(h(i,j) <= 1.)fsm(i,j)=0.

            if(fsm(i,j) == 1.)totarea=totarea+dx(i,j)*dy(i,j)*1.e-6
        enddo
    enddo
!---find coast
    totcoast=0.
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
                    totcoast=totcoast+10000.
                endif
            endif
        enddo
    enddo
    write(*,*)totcoast
!    stop
    write(*,*)'ZZ',zz(kb-1)
    directory='../data/HindPlastic31/FISHOUT/'
    icount=0
!_--open
!    open(133,file='MacroValida_2010_2012_2.dat')
!    open(144,file='MicroValida_2010_2012_2.dat')
    open(123,file='ibm2010_2012_bur_x05_retx2.dat')
!    open(144,file='popdepav_microSize.dat')
!    open(11,file='popAll2010seasDep1.bin',form='unformatted')
    open(12,file='popAll2010seasDep10.bin',form='unformatted') ! auto thelw!!!
!    open(13,file='depAll2010seas.bin',form='unformatted')
!    open(15,file='popAll2010seasBot.bin',form='unformatted')
!    open(20,file='popAll2010seasDep10biom.bin'form='unformatted')
!    open(19,file='daysAll2010seasDep10.bin',form='unformatted')
!    open(21,file='distAll2010seasDep10.bin',form='unformatted')
!    open(14,file='popAll2010seasBea.bin',form='unformatted')
!    open(21,file='pop3dAll2010seas.bin',form='unformatted')
!    open(16,file='denAll2010seasDep10.bin',form='unformatted')
!    open(17,file='sinAll2010seasDep10.bin',form='unformatted')
!    open(18,file='fouAll2010seasDep10.bin',form='unformatted')

    unit1=1.E-6
    ldf=0
    do ii=1,120
        chara=directory(ii:ii)
        if(chara /= ' ')then
            ldf=ldf+1
        endif
    enddo

    fname=directory(1:ldf)//'listdays'
    ldf1=ldf+13
    write(*,*)fname(1:ldf1)

    open(33,file=fname(1:ldf1),status='old')
!    do iday=1,730
!        read(33,'(a26)',end=3000)fname1
!    enddo

    ndays=3  !31   !365*3
!    ndays=10
    nsum=3  !31   !365*3 !edw vazw tis meres pou exei treksei
!    nsum=1
    do 20 iday=1,ndays
        read(33,'(a26)',end=3000)fname1
        write(*,*)fname1
        fname2=directory(1:ldf)//fname1
        ldf2=ldf+20
        write(*,*)fname2(1:ldf2)
        iwrite=0
        if(mod(iday,nsum) == 0)iwrite=1
        goto 789
        status = NF_OPEN(fname2(1:ldf2),0, ncid)
        write(*,*)'OK',filename,ncid
    !---get variables
        status=nf_inq_varid(ncid,'species',species_id)
        status=nf_inq_varid(ncid,'nfi',nfi_id)
        status=nf_inq_varid(ncid,'variant',variant_id)
        status=nf_inq_varid(ncid,'xsi',xsi_id)
        status=nf_inq_varid(ncid,'ysi',ysi_id)
        status=nf_inq_varid(ncid,'ref',ref_id)
        status=nf_inq_varid(ncid,'iarea',iarea_id)
        status=nf_inq_varid(ncid,'fdep',fdep_id)
        status=nf_inq_varid(ncid,'rbt',rbt_id)
        status=nf_inq_varid(ncid,'fou',fou_id)
    !    status=nf_inq_varid(ncid,'icountdays',icountdays_id)
        status=nf_get_var_int(ncid,species_id,species)
        status=nf_get_var_int(ncid,variant_id,variant)
        status=nf_get_var_int(ncid,iarea_id,iarea)
        status=nf_get_var_int(ncid,ref_id,ref)
    !    status=nf_get_var_int(ncid,icountdays_id,icountdays)
        status=nf_get_var_real(ncid,nfi_id,nfi)
        status=nf_get_var_real(ncid,xsi_id,xsi)
        status=nf_get_var_real(ncid,ysi_id,ysi)
        status=nf_get_var_real(ncid,rbt_id,rbt)
        status=nf_get_var_real(ncid,fou_id,fou)
        status=nf_get_var_real(ncid,fdep_id,fdep)
        status = NF_CLOSE (NCID)
        789 continue
        write(*,*)'NSI',nsi
        nsi=nmax
        open(99,file=fname2(1:ldf2),form='unformatted',status='old')
        read(99)nsi
        read(99)(species(n),n=1,nsi)
    !    read(99)(xfi(n),n=1,nsi)
        read(99)(nfi(n),n=1,nsi)
        read(99)(variant(n),n=1,nsi)
        read(99)(xsi(n),n=1,nsi)
        read(99)(ysi(n),n=1,nsi)
        read(99)(ref(n),n=1,nsi)
    !    read(99)(iarea(n),n=1,nsi)
        read(99)(fdep(n),n=1,nsi)
        read(99)(icountdays(n),n=1,nsi)
        read(99)(fou(n),n=1,nsi)
        read(99)(dista(n),n=1,nsi)
        close(99)

        nsi=nmax
        write(*,*)'NSI',nsi
        do k=1,nsi
            if(species(k) == 1. .AND. variant(k) > 4.)then
                fdep(k)=0.
            endif
        enddo

        icount=icount+1
        concbeach3=0.
        icountbotbea=0
    !    nsum=nsum+1
        ic=0
        do 30 ch=1,nsi
            factor=1.
    !        if(iarea(ch).eq.1.)then
    !            if(species(ch).eq.2.and.variant(ch).eq.2)factor=0.25
    !            if(species(ch).eq.2.and.variant(ch).eq.1)factor=1.25
    !        endif
            if(species(ch) == 1 .AND. variant(ch) == 1)factWei=1.25e-7
            if(species(ch) == 1 .AND. variant(ch) == 2)factWei=1.e-3
            if(species(ch) == 2 .AND. variant(ch) == 1)factWei=0.1
            if(species(ch) == 2 .AND. variant(ch) == 2)factWei=1.
            if(species(ch) == 2 .AND. variant(ch) == 3)factWei=15.
            if(species(ch) == 2 .AND. variant(ch) == 4)factWei=5.
            if(species(ch) == 2 .AND. variant(ch) == 5)factWei=4.

            rx= rand(0)
            rr=(2. * rx - 1.)
            ii=nint((xsi(ch)-xleft)/xs)+1
            jj=nint((ysi(ch)-ybot)/ys)+1

            if(variant(ch) == .0 .OR. species(ch) == 0)goto 123

            ic=ic+1
            isiind(ic)=ch
            nsp=species(ch)
            ivariant=variant(ch)
            iref=ref(ch)
               
            if(ii > im .OR. ii < 1 .OR. jj > jm .OR. jj < 1)goto 123

        !---biomass per type (tn/Km2?)
            concbiom(ii,jj,ivariant,nsp)=concbiom(ii,jj,ivariant,nsp)+ &
            factor*nfi(ch)*xfiInit(ivariant,nsp)*1.e-6/unitArea/float(nsum)

        !---total biomass
            if(fdep(ch) > -11. .AND. iref == 0)then
                concbt(ii,jj)=concbt(ii,jj)+ &
                factor*nfi(ch)/unitArea/float(nsum)*factWei
            endif
                  
        !---concentration (#/Km2)
            if(fdep(ch) > -11. .AND. iref == 0)then
                conc(ii,jj,ivariant,nsp)=conc(ii,jj,ivariant,nsp)+ &
                factor*nfi(ch)/unitArea/float(nsum)
            endif
            if(fdep(ch) < -600. .AND. iref == 0)then
                conc2(ii,jj,ivariant,nsp)=conc2(ii,jj,ivariant,nsp)+ &
                factor*nfi(ch)/unitArea/float(nsum)
            endif
            if(fdep(ch) > -0.3 .AND. iref == 0)then
                conc1(ii,jj,ivariant,nsp)=conc1(ii,jj,ivariant,nsp)+ &
                factor*nfi(ch)/unitArea/float(nsum)
            endif
            if(nsp == 1 .OR. (nsp == 2 .AND. ivariant == 4))then
            !---3d conc
                if(nsp == 1 .AND. iref == 0)then
                    do k=2,klev2
                        ddpt=dpt2(k)-dpt2(k-1)
                        if(-fdep(ch) >= dpt2(k-1) .AND. -fdep(ch) <= dpt2(k)) then
                        !---convert to #/m3
                            conc3d(ii,jj,k-1,ivariant)=conc3d(ii,jj,k-1,ivariant)+ &
                            factor*nfi(ch)/unitArea/float(nsum)/ddpt*1.e-6
                        endif
                    enddo
                endif
            !---calculate sinking velocity and density
                rhom = 1025.
                den = dens(ivariant,nsp) !particle density
                denf=1500.
                dia0 = diam(ivariant,nsp) !particle diameter
                dia02= 0.5*dia0       !particle radius
                den0 = dens(ivariant,nsp) !particle density
                hfou=fou(ch)
                dia2 = dia02+hfou
                dia  = 2.*dia2
                if(nsp == 1)then
                    hth=dia02*(((denf-den0)/(denf-rhom))**0.333-1)
                else
                    hth=thic*(rhom-den0)/(denf-rhom)*0.5
                endif
                if(nsp == 2 .AND. ivariant == 4)then
                    den =den0*thic/(thic+2*hfou)+ &
                    denf*2*hfou/(thic+2*hfou)
                    cd=0.3
                else
                    den=(den0*dia02**3.+denf*(dia2**3.-dia02**3.))/(dia2**3.)
                    cd=1.
                endif
                densa(ch)=den
            !---buoyancy
                grav=9.806e0
                dc = 9.52 * vk ** (2. / 3.) / &
                (grav** (1. /3.) * (1. - den /rhom) ** (1. / 3.))
                if(dia <= dc) then
                    wb = cd*grav * dia ** 2 * (1. - den / rhom) / (18. * vk)
                else
                    wb = cd*8. / 3. * grav * dia * sign(1., 1. - den / rhom) &
                    * sqrt(abs(1. - den / rhom))
                end if
                if(nsp == 1 .AND. ivariant == 1 .AND. ref(ch) == 0)then
                    id=icountdays(ch)
                    fdepav1(id)=fdepav1(id)+fdep(ch)*nfi(ch)
                    wbav1(id)=wbav1(id)+wb*nfi(ch)
                    hfouav1(id)=hfouav1(id)+fou(ch)*nfi(ch)
                    denav1(id)=denav1(id)+den*nfi(ch)
                    count1(id)=count1(id)+nfi(ch)
                endif
                if(nsp == 1 .AND. ivariant == 2 .AND. ref(ch) == 0)then
                    id=icountdays(ch)
                    fdepav2(id)=fdepav2(id)+fdep(ch)*nfi(ch)
                    wbav2(id)=wbav2(id)+wb*nfi(ch)
                    hfouav2(id)=hfouav2(id)+fou(ch)*nfi(ch)
                    denav2(id)=denav2(id)+den*nfi(ch)
                    count2(id)=count2(id)+nfi(ch)
                endif
                if(fdep(ch) > -500.)then
                    sink2d(ii,jj,ivariant,nsp)=sink2d(ii,jj,ivariant,nsp)+wb
                    density(ii,jj,ivariant,nsp)=density(ii,jj,ivariant,nsp)+den
                    countSink(ii,jj,ivariant,nsp)=countSink(ii,jj,ivariant,nsp)+1
                endif
            endif
        !---end sinking
                   
        !---mean depth
            sidep(ii,jj,ivariant,nsp)=sidep(ii,jj,ivariant,nsp)+fdep(ch)*nfi(ch)
            countsi1(ii,jj,ivariant,nsp)=countsi1(ii,jj,ivariant,nsp)+nfi(ch)
            if(fdep(ch) > -10.)then
            !---mean days at sea
                sidays(ii,jj,ivariant,nsp)=sidays(ii,jj,ivariant,nsp)+ &
                icountdays(ch)*nfi(ch)
                if(dista(ch) < 1.e+6)then
                    sidist(ii,jj,ivariant,nsp)=sidist(ii,jj,ivariant,nsp)+ &
                    dista(ch)*nfi(ch)
                endif
                countsi2(ii,jj,ivariant,nsp)=countsi2(ii,jj,ivariant,nsp)+nfi(ch)
            !---mean biofilm
                sifou(ii,jj,ivariant,nsp)=sifou(ii,jj,ivariant,nsp)+fou(ch)/hth*nfi(ch)
            endif

        !---beach conc.
            if(ref(ch) == 77 .AND. coast(ii,jj) == 1)then
                concbeach(ii,jj,ivariant,nsp)=concbeach(ii,jj,ivariant,nsp)+ &
                factor*nfi(ch)/unitArea/float(nsum)
                if(ivariant == 3 .AND. nsp == 2)then
                    concbeach3=concbeach3+nfi(ch)/5000.
                    icountbotbea=icountbotbea+1
                endif
            endif

        !---bottom conc.
            if(ref(ch) == -1 .OR. (ref(ch) == 0 .AND. &
            fdep(ch) < -0.9*h(ii,jj)))then
                concbottom(ii,jj,ivariant,nsp)=concbottom(ii,jj,ivariant,nsp)+ &
                factor*nfi(ch)/unitArea/float(nsum)
            endif
        123 continue
        30 END DO

        do ivariant=1,6
            do k=2,klev2
                icount=0
                do i=1,im
                    do j=1,jm
                        if(conc3d(i,j,k-1,ivariant) > 0.) then
                            conc3dav(k-1,ivariant)=conc3dav(k-1,ivariant)+conc3d(i,j,k-1,ivariant)
                            icount=icount+1
                        endif
                    enddo
                enddo
                if(icount > 0)then
                    conc3dav(k-1,ivariant)=conc3dav(k-1,ivariant)/float(icount)
                endif
            enddo
        enddo

        reftot=0
        beach=0
        beach2=0
        siback=0
        sisource=0
        sibottom=0
        empty=0
        do ch=1,nsi
            ii=nint((xsi(ch)-xleft)/xs)+1
            jj=nint((ysi(ch)-ybot)/ys)+1
             
        !---calculate total gone out of model area
            if(variant(ch) == 0)empty=empty+1
            if(variant(ch) > 0 .AND. species(ch) > 0)then
                if(ref(ch) == -2)then
                    reftot=reftot+1
                endif
                if(ref(ch) == 77)then
                    beach=beach+1
                endif
                if(ref(ch) == 99)then
                    beach2=beach2+1
                endif
                if(iarea(ch))then
                    siback=siback+1
                endif
                if(iarea(ch) > 1)then
                    sisource=sisource+1
                endif
                if(ref(ch) == -1)then
                    sibottom=sibottom+1
                endif
            endif
        !----------------------------------------
        enddo ! end of Sis

        write(*,*)'SIs out-of-domain',reftot
        write(*,*)'SIs background',siback
        write(*,*)'SIs sources',sisource
        write(*,*)'SIs beach',beach
        write(*,*)'SIs beach2',beach2
        write(*,*)'SIs empty',empty
        write(*,*)'SIs bottom',sibottom
        call totfish(1,totarea,totcoast)

        do nsp=1,nspecies
            totbiom(nsp)=0.
            totbiomBea(nsp)=0.
            totbiom100(nsp)=0.
            totbiomBot(nsp)=0.
            do ivariant=1,nvariantsmax
                totbiom100(nsp)=totbiom100(nsp)+tbvariant100(ivariant,nsp)*1.e-6
                totbiom(nsp)=totbiom(nsp)+tbvariant(ivariant,nsp)*1.e-6
                totbiomBea(nsp)=totbiomBea(nsp)+ &
                (tbvariantBea(ivariant,nsp)+tbvariantBur(ivariant,nsp))*1.e-6
                totbiomBot(nsp)=totbiomBot(nsp)+ &
                tbvariantBot(ivariant,nsp)*1.e-6
                concmean(ivariant,nsp)=concmean(ivariant,nsp)+ &
                tnvariant100(ivariant,nsp)
                concmeanBea(ivariant,nsp)=concmeanBea(ivariant,nsp)+ &
                tnvariantBea(ivariant,nsp)
            enddo
        enddo

        write(123,7777) &
        (tnvariant100(ivariant,1),ivariant=1,nvariants(1)), &
        (tnvariant100(ivariant,2),ivariant=1,nvariants(2)), &
        (tnvariantBea(ivariant,2),ivariant=1,nvariants(2)), &
        (tbvariant(ivariant,1)*1.e-6,ivariant=1,nvariants(1)), &
        (tbvariant(ivariant,2)*1.e-6,ivariant=1,nvariants(2)), &
        (tbvariant100(ivariant,1)*1.e-6,ivariant=1,nvariants(1)), &
        (tbvariant100(ivariant,2)*1.e-6,ivariant=1,nvariants(2)), &
        ((tbvariantBea(ivariant,1)+tbvariantBur(ivariant,1))*1.e-6,ivariant=1,nvariants(1)), &
        ((tbvariantBea(ivariant,2)+tbvariantBur(ivariant,2))*1.e-6,ivariant=1,nvariants(2)), &
        (tbvariantBot(ivariant,1)*1.e-6,ivariant=1,nvariants(1)), &
        (tbvariantBot(ivariant,2)*1.e-6,ivariant=1,nvariants(2)), &
        (fdepvariant(ivariant,1),ivariant=1,nvariants(1)), &
        fdepvariant(4,2),fdepvariant(3,2), &
        (fdepvariantBot(ivariant,1),ivariant=1,nvariants(1)), &
        fdepvariantBot(4,2),fdepvariantBot(3,2)   !  &
    !     (denvariant(ivariant,1),ivariant=1,nvariants(1)), &
    !     denvariant(4,2),denvariant(3,2)

        write(*,*)'BIOMASS-MICRO',totbiom(1),totbiom100(1),totbiomBea(1), &
        totbiomBot(1)
        write(*,*)'BIOMASS-MACRO',totbiom(2),totbiom100(2),totbiomBea(2), &
        totbiomBot(2)

!       write(123,7777)(tbvariant100(2,1)+tbvariant100(1,1))*1.e-6, &
!      (tbvariant100(1,2)+tbvariant100(2,2)+tbvariant100(3,2)+          &
!       tbvariant100(4,2)+tbvariant100(5,2))*1.e-6,                 &
!           tnvariant100(1,1),tnvariant100(2,1),tnvariant100(1,2),      &
!           tnvariant100(2,2),tnvariant100(3,2),tnvariant100(4,2),      &
!           tnvariant100(5,2),                                  &
!           1.e-6*tbvariant100(1,1),1.e-6*tbvariant100(2,1),        &
!           1.e-6*tbvariant100(1,2),1.e-6*tbvariant100(2,2),        &
!       1.e-6*tbvariant100(3,2),1.e-6*tbvariant100(4,2),1.e-6*tbvariant100(5,2), &
!       1.e-6*tbvariantBea(1,1),1.e-6*tbvariantBea(2,1),            &
!           1.e-6*tbvariantBea(1,2),1.e-6*tbvariantBea(2,2),        &
!       1.e-6*tbvariantBea(3,2),1.e-6*tbvariantBea(4,2),1.e-6*tbvariantBea(5,2), &
!       1.e-6*tbvariantBot(1,1),1.e-6*tbvariantBot(2,1),1.e-6*tbvariantBot(4,2), 

        7777 format(105(e12.5,2X))

        if(iwrite == 1)then
            write(*,*)'iwrite'
            call totfish(1,totarea,totcoast)
            write(*,*)'totfish'
            do nsp=1,nspecies
                do ivariant=1,nvariantsmax
                    do i=1,im2
                        do j=1,jm2
                            sidep(i,j,ivariant,nsp)=sidep(i,j,ivariant,nsp)/countsi1(i,j,ivariant,nsp)
                            sidays(i,j,ivariant,nsp)=sidays(i,j,ivariant,nsp)/countsi2(i,j,ivariant,nsp)
                            sifou(i,j,ivariant,nsp)=sifou(i,j,ivariant,nsp)/countsi2(i,j,ivariant,nsp)
                            sidist(i,j,ivariant,nsp)=sidist(i,j,ivariant,nsp)/countsi2(i,j,ivariant,nsp)
                        enddo
                    enddo
                enddo
            enddo
            do nsp=1,nspecies
                do ivariant=1,nvariantsmax
                    do i=1,im2
                        do j=1,jm2
                            sink2d(i,j,ivariant,nsp)=sink2d(i,j,ivariant,nsp)/ &
                            countsink(i,j,ivariant,nsp)
                            density(i,j,ivariant,nsp)=density(i,j,ivariant,nsp)/ &
                            countsink(i,j,ivariant,nsp)
                        enddo
                    enddo
                enddo
            enddo
        !---write monthly mean fields for eggs/adults etc
            do nsp=1,nspecies
                do ivariant=1,nvariants(nsp)
                    if(nsp == 1)then
                        do k=1,klev2-1
                            if(conc3dav(k,ivariant) == 0.)conc3dav(k,ivariant)=-1.E+10
                        enddo
                    endif
                    do j=1,jm
                        do i=1,im
                    !       if(nsp.eq.1.and.ivariant.eq.5) &
                    !           write(*,*)'TEST1',i,j,conc(i,j,ivariant,nsp)
                            if(conc(i,j,ivariant,nsp) == 0.)conc(i,j,ivariant,nsp)=-1.E+10
                            if(conc1(i,j,ivariant,nsp) == 0.)conc1(i,j,ivariant,nsp)=-1.E+10
                            if(conc2(i,j,ivariant,nsp) == 0.)conc2(i,j,ivariant,nsp)=-1.E+10
                            if(sidep(i,j,ivariant,nsp) == 0.)sidep(i,j,ivariant,nsp)=-1.E+10
                            if(concbeach(i,j,ivariant,nsp) == 0.)concbeach(i,j,ivariant,nsp)=-1.E+10
                            if(concbottom(i,j,ivariant,nsp) == 0.)concbottom(i,j,ivariant,nsp)=-1.E+10
                            if(sidays(i,j,ivariant,nsp) == 0.)sidays(i,j,ivariant,nsp)=-1.E+10
                            if(sidist(i,j,ivariant,nsp) == 0.)sidist(i,j,ivariant,nsp)=-1.E+10
                            if(sifou(i,j,ivariant,nsp) == 0.)sifou(i,j,ivariant,nsp)=-1.E+10
                            if(nsp == 1)then
                                do k=1,klev2-1
                                    if(conc3d(i,j,k,ivariant) == 0.)conc3d(i,j,k,ivariant)=-1.E+10
                                enddo
                            endif
                        enddo
                    enddo
                    do j=1,jm
                        do i=1,im
                            if(concbt(i,j) == 0.)concbt(i,j)=-1.E+10
                        enddo
                    enddo
                    write(10)(conc2(i,1,ivariant,nsp),i=1,im2*jm2)
                    write(12)(conc(i,1,ivariant,nsp),i=1,im2*jm2)
                !    if(nsp.eq.1.and.ivariant.ge.3)then
                    write(11)(conc1(i,1,ivariant,nsp),i=1,im2*jm2)
                !    endif
                    write(13)(sidep(i,1,ivariant,nsp),i=1,im2*jm2)
                    write(14)(concbeach(i,1,ivariant,nsp),i=1,im2*jm2)
                    write(15)(concbottom(i,1,ivariant,nsp),i=1,im2*jm2)
                !    write(16)(density(i,1,ivariant,nsp),i=1,im2*jm2)
                    write(17)(sink2d(i,1,ivariant,nsp),i=1,im2*jm2)
                    write(18)(sifou(i,1,ivariant,nsp),i=1,im2*jm2)
                    write(19)(sidays(i,1,ivariant,nsp),i=1,im2*jm2)
                    write(21)(sidist(i,1,ivariant,nsp),i=1,im2*jm2)
                enddo
            enddo

            icountbotbea=0
            concbeach3=0.
            do i=1,im
                do j=1,jm
                    if(concbeach(i,j,3,2) > 0.)then
                        concbeach3=concbeach3+concbeach(i,j,3,2)
                        icountbotbea=icountbotbea+1
                    endif
                enddo
            enddo
            concbeach3=concbeach3/float(icountbotbea)*5/1000.
            write(*,*)'BOTBEAM',concbeach3
    !        write(20)(concbt(i,1),i=1,im2*jm2)
    !        write(21)(conc3d(i,1,1,5),i=1,im2*jm2)
    !        write(21)(conc3d(i,1,1,6),i=1,im2*jm2)
                   
            do ivariant=1,4
                do k=1,klev2-1
            !       write(21)(conc3d(i,1,k,ivariant),i=1,im2*jm2)
            !       write(12)(conc3d(i,1,k,2),i=1,im2*jm2)
                enddo
            enddo

        !---Model validation
            open(1,file='macroAllnew.dat')
            open(2,file='microAllnew300_2.dat')
            write(*,*)'VALIDATION MACRO',nobsmacro
            do nn=1,nobsmacro
                read(1,*)xlon,xlat
                write(*,*)nn,xlat,xlon
                ii=nint((xlon-xleft)/xs)+1
                jj=nint((xlat-ybot)/xs)+1
                iseas=(imo-1)/3+1
                do ivariant=1,nvariants(2)
                    if(conc(ii,jj,ivariant,2) <= 0.)then
                        val=0.
                        icount=0
                        do i=ii-2,ii+2
                            do j=jj-2,jj+2
                                if(conc(i,j,ivariant,2) > 0.)then
                                    val=val+conc(i,j,ivariant,2)
                                    icount=icount+1
                                endif
                            enddo
                        enddo
                        if(icount > 0)then
                            val=val/float(icount)
                            conc_model_macro(ivariant)=val
                        else
                            conc_model_macro(ivariant)=0.
                        endif
                    else
                        conc_model_macro(ivariant)=conc(ii,jj,ivariant,2)
                    endif
                enddo
                write(133,*)iseas,conc_model_macro(1),conc_model_macro(2), &
                conc_model_macro(3)+conc_model_macro(4)+conc_model_macro(5)
            enddo

            write(*,*)'VALIDATION MICRO',nobsmicro

            do n=1,nobsmicro
                read(2,*)iyear,imo,iiday,xlon,xlat
                write(*,*)n,xlat,xlon
                ii=nint((xlon-xleft)/xs)+1
                jj=nint((xlat-ybot)/xs)+1
                iseas=(imo-1)/3+1
                do ivariant=1,nvariants(1)
                    if(conc1(ii,jj,ivariant,1) <= 0.)then
                        val=0.
                        icount=0
                        do i=ii-2,ii+2
                            do j=jj-2,jj+2
                                if(conc1(i,j,ivariant,1) > 0.)then
                                    val=val+conc1(i,j,ivariant,1)
                                    icount=icount+1
                                endif
                            enddo
                        enddo
                        if(icount > 0)then
                            val=val/float(icount)
                            conc_model_micro(ivariant)=val
                        else
                            conc_model_micro(ivariant)=0.
                        endif
                    else
                        conc_model_micro(ivariant)=conc1(ii,jj,ivariant,1)
                    endif
                enddo
                write(144,*)iseas,conc_model_micro(3),conc_model_micro(4), &
                conc_model_micro(5),conc_model_micro(6)
            enddo
        !---end validation

            do nsp=1,nspecies
                if(nsp == 1)then
                    do k=1,klev2-1
                !        write(144,189)(conc3dav(k,ivariant),ivariant=1,nvariants(nsp))
                        189 format(10(e12.5,2X))
                    enddo
                endif
                do ivariant=1,nvariants(nsp)
                    if(nsp == 1)then
                        do k=1,klev2-1
                            conc3dav(k,ivariant)=0.
                        enddo
                    endif
                    do i=1,im2
                        do j=1,jm2
                            conc1(i,j,ivariant,nsp)=0.
                            conc2(i,j,ivariant,nsp)=0.
                            conc(i,j,ivariant,nsp)=0.
                            concbt(i,j)=0.
                            concbiom(i,j,ivariant,nsp)=0.
                            concbeach(i,j,ivariant,nsp)=0.
                            concbottom(i,j,ivariant,nsp)=0.
                            sidep(i,j,ivariant,nsp)=0.
                            countsi2(i,j,ivariant,nsp)=0.
                            countsi1(i,j,ivariant,nsp)=0.
                            countsink(i,j,ivariant,nsp)=0.
                            sink2d(i,j,ivariant,nsp)=0.
                            density(i,j,ivariant,nsp)=0.
                            sidays(i,j,ivariant,nsp)=0.
                            sidist(i,j,ivariant,nsp)=0.
                            if(nsp == 1)then
                                do k=1,klev2-1
                                    conc3d(i,j,k,ivariant)=0.
                                enddo
                            endif
                        enddo
                    enddo
                enddo
            enddo

    !        write(*,*)'NSUM',nsum,conc(130,130,1,1)
    !        nsum=0


    !        do id=1,365
    !            if(count1(id).gt.0)then
    !                fdepav1(id)=fdepav1(id)/count1(id)
    !                wbav1(id)=wbav1(id)/count1(id)
    !                 hfouav1(id)=hfouav1(id)/count1(id)
    !                 denav1(id)=denav1(id)/count1(id)
    !                 write(143,*)id,fdepav1(id),wbav1(id),hfouav1(id),denav1(id)
    !             else
    !                 write(143,*)id,'NaN ','NaN ','NaN ','NaN '
    !             endif
    !             if(count2(id).gt.0)then
    !                 fdepav2(id)=fdepav2(id)/count2(id)
    !                 wbav2(id)=wbav2(id)/count2(id)
    !                 hfouav2(id)=hfouav2(id)/count2(id)
    !                 denav2(id)=denav2(id)/count2(id)
    !                 write(144,*)id,fdepav2(id),wbav2(id),hfouav2(id),denav2(id)
    !             else
    !                 write(144,*)id,'NaN ','NaN ','NaN ','NaN '
    !             endif
    !         enddo
            
    !         do id=1,365
    !             fdepav1(id)=0.
    !             wbav1(id)=0.
    !             hfouav1(id)=0.
    !             denav1(id)=0.
    !             count1(id)=0.
    !             fdepav2(id)=0.
    !             wbav2(id)=0.
    !             hfouav2(id)=0.
    !             denav2(id)=0.
    !             count2(id)=0.
    !         enddo
        endif
        concbeach3=concbeach3/float(icountbotbea)
        write(*,*)'BOTBEA',concbeach3
    20 enddo
    3000 continue
    do nsp=1,nspecies
        do ivariant=1,nvariants(nsp)
            concmeanBea(ivariant,nsp)=concmeanBea(ivariant,nsp)/float(icount)
            concmean(ivariant,nsp)=concmean(ivariant,nsp)/float(icount)
            write(*,*)'MEAN CONC',nsp,ivariant,concmean(ivariant,nsp), &
            concmeanBea(ivariant,nsp)
        enddo
    enddo
    stop
end program





subroutine totfish(iwrite,totarea,totcoast)
!! ??
    include 'ibmall.fcm'
    integer :: ch,i,nn,iwres
    integer :: icountempty,iwrite
    real :: weira(nvariantsmax,nspecies)
    integer :: countdays(nvariantsmax,nspecies)

    do nsp=1,nspecies
        do i=1,nvariants(nsp)
            weivariant(i,nsp)=0.
            tbvariant(i,nsp)=0.
            tbvariant100(i,nsp)=0.
            tbvariantBea(i,nsp)=0.
            tbvariantBot(i,nsp)=0.
            tnvariantBea(i,nsp)=0.
            tnvariantBot(i,nsp)=0.
            tnvariant(i,nsp)=0.
            tnvariant100(i,nsp)=0.
            tsivariant(i,nsp)=0.
            tsivariant100(i,nsp)=0.
            tsivariantBot(i,nsp)=0.
            tsivariantBea(i,nsp)=0.
            denvariant(i,nsp)=0.
            fdepvariant(i,nsp)=0.
            fdepvariantBot(i,nsp)=0.
            countdays(i,nsp)=0.
        enddo
    enddo
!---Calculation of empty SIs,then fill with new eggs in ersem_fishGrad.F
    do n=1,nmax
        nempty(n)=0
    enddo
    icountempty=0
    icountsource=0
    icountback=0
    do n=1,nsi
        factor=1.
!        if(species(n).eq.2.and.variant(n).ge.3)factor=2.
!        if(species(n).eq.2.and.variant(n).lt.3)factor=5.
        if(variant(n) == 0)then
            icountempty=icountempty+1
            nempty(icountempty)=n
        else
            if(iarea(n) > 1)icountsource=icountsource+1
            if(iarea(n) == 1)icountback=icountback+1
            nsp=species(n)
            ivariant=variant(n)
            if(ref(n) == 77)then
                tnvariantBea(ivariant,nsp)=tnvariantBea(ivariant,nsp)+nfi(n)
                if(nsp == 1)factBea=0.075
                if(nsp == 2 .AND. ivariant <= 2)factBea=0.02
                if(nsp == 2 .AND. ivariant > 2)factBea=0.05
                tbvariantBea(ivariant,nsp)=tbvariantBea(ivariant,nsp)+ &
                xfiInit(ivariant,nsp)*nfi(n)*factor
                tbvariantBur(ivariant,nsp)=tbvariantBur(ivariant,nsp)+ &
                xfiInit(ivariant,nsp)*nfi(n)*factBea
                tsivariantBea(ivariant,nsp)=tsivariantBea(ivariant,nsp)+1
            endif
            if(ref(n) == -1)then
                tnvariantBot(ivariant,nsp)=tnvariantBot(ivariant,nsp)+ &
                nfi(n)*factor
                tbvariantBot(ivariant,nsp)=tbvariantBot(ivariant,nsp)+ &
                xfiInit(ivariant,nsp)*nfi(n)*factor
                tsivariantBot(ivariant,nsp)=tsivariantBot(ivariant,nsp)+1
                fdepvariantBot(ivariant,nsp)=fdepvariantBot(ivariant,nsp)+fdep(n)*nfi(n)
            endif
            if(ref(n) == 0)then
                tbvariant(ivariant,nsp)=tbvariant(ivariant,nsp)+ &
                nfi(n)*xfiInit(ivariant,nsp)*factor
                if(fdep(n) > -10)then
                    tbvariant100(ivariant,nsp)=tbvariant100(ivariant,nsp)+ &
                    nfi(n)*xfiInit(ivariant,nsp)*factor
                    tnvariant100(ivariant,nsp)=tnvariant100(ivariant,nsp)+ &
                    nfi(n)*factor
                    tsivariant100(ivariant,nsp)=tsivariant100(ivariant,nsp)+1
                endif
                if(isnan(densa(n)) == .FALSE. .AND. &
                isnan(fdep(n)) == .FALSE. .AND. ref(n) == 0)then
                    tnvariant(ivariant,nsp)=tnvariant(ivariant,nsp)+ &
                    nfi(n)*factor
                endif
                tsivariant(ivariant,nsp)=tsivariant(ivariant,nsp)+1

                weivariant(ivariant,nsp)=weivariant(ivariant,nsp)+ &
                nfi(n)*xfiInit(ivariant,nsp)
!                if(isnan(fdep(n)))write(*,*)'TESTBAG',fou(n),ref(n),nfi(n)
                if(isnan(fdep(n)) == .FALSE. .AND. ref(n) == 0)then
                    fdepvariant(ivariant,nsp)=fdepvariant(ivariant,nsp)+fdep(n)*nfi(n)
                    denvariant(ivariant,nsp)=denvariant(ivariant,nsp)+densa(n)*nfi(n)
                endif
                countdays(ivariant,nsp)=countdays(ivariant,nsp)+ &
                icountdays(n)*nfi(n)
            endif
        endif
    enddo

    do n=1,nspecies
        do i=1,nvariants(n)
            if(tnvariant(i,n) > 0.)then
                weivariant(i,n)=weivariant(i,n)/tnvariant(i,n)
                denvariant(i,n)=denvariant(i,n)/tnvariant(i,n)
                fdepvariant(i,n)=fdepvariant(i,n)/tnvariant(i,n)
                countdays(i,n)=countdays(i,n)/tnvariant(i,n)
            else
                weivariant(i,n)=0.
                fdepvariant(i,n)=0.
                countdays(i,n)=0.
                denvariant(i,n)=0.
            endif
            if(tnvariantBot(i,n) > 0.)then
                fdepvariantBot(i,n)=fdepvariantBot(i,n)/tnvariantBot(i,n)
            else
                fdepvariantBot(i,n)=0.
            endif
            tnvariant100(i,n)=tnvariant100(i,n)/totarea
            tbvariant100(i,n)=tbvariant100(i,n)
            tnvariant(i,n)=tnvariant(i,n)/totarea
!            tnvariantBea(i,n)=tnvariantBea(i,n)/tsivariantBea(i,n)
            tnvariantBea(i,n)=tnvariantBea(i,n)/totcoast
            tbvariantBea(i,n)=tbvariantBea(i,n)
            tnvariantBot(i,n)=tnvariantBot(i,n)/totarea
            tbvariantBot(i,n)=tbvariantBot(i,n)
        enddo
    enddo

    if(iwrite == 1)then
!        write(*,*)'NSI=',nsi
!        write(*,*)'EMPTY=',icountempty
!        write(*,*)'SIs Source',icountsource
!        write(*,*)'SIs Background',icountback
!        do n=1,nspecies
!            write(*,*)'SPECIES',n
        adultbiom=0.
!            do ivariant=1,nvariants(n)
!                write(*,*)'TOTAL BIOMASS',ivariant,tbvariant(ivariant,n)*1.e-6
!                write(*,*)'PARTICLES',ivariant,tnvariant(ivariant,n)
!                write(*,*)'SIs',ivariant,tsivariant(ivariant,n)
!                write(*,*)'PARTICLE WEIGHT(mgr)',ivariant,weivariant(ivariant,n)
!                write(*,*)'DEPTH(m)',ivariant,fdepvariant(ivariant,n)
!                write(*,*)'DAYS-AT-SEA',ivariant,countdays(ivariant,n)
!                write(*,*)'DENS',ivariant,denvariant(ivariant,n)
!            enddo
!        enddo
        write(*,*)'SIs ',' SURF-CONC ',' BEA-CONC ',' BOT-CONC '
        write(*,*)'MICRO1',tsivariant(1,1), &
        tnvariant100(1,1),tnvariantBea(1,1),tbvariantBot(1,1)
        write(*,*)'MICRO2',tsivariant(2,1), &
        tnvariant100(2,1),tnvariantBea(2,1),tbvariantBot(2,1)
        write(*,*)'MICRO3',tsivariant(3,1), &
        tnvariant100(3,1),tnvariantBea(3,1),tbvariantBot(3,1)
        write(*,*)'MICRO4',tsivariant(4,1), &
        tnvariant100(4,1),tnvariantBea(4,1),tbvariantBot(4,1)
        write(*,*)'MICRO5',tsivariant(5,1), &
        tnvariant100(5,1),tnvariantBea(5,1),tbvariantBot(5,1)
        write(*,*)'MICRO6',tsivariant(6,1), &
        tnvariant100(6,1),tnvariantBea(6,1),tbvariantBot(6,1)
!        write(*,*)'MACRO1',tsivariant(1,2),
        write(*,*)'MACRO1',tnvariant(1,2), &
        tnvariant100(1,2),tnvariantBea(1,2),tbvariantBot(1,2)
        write(*,*)'MACRO2',tsivariant(2,2), &
        tnvariant100(2,2),tnvariantBea(2,2),tbvariantBot(2,2)
        write(*,*)'MACRO3',tsivariant(3,2), &
        tnvariant100(3,2),tnvariantBea(3,2),tbvariantBot(3,2)
        write(*,*)'MACRO4',tsivariant(4,2), &
        tnvariant100(4,2),tnvariantBea(4,2),tbvariantBot(4,2)
        write(*,*)'MACRO5',tsivariant(5,2), &
        tnvariant100(5,2),tnvariantBea(5,2),tbvariantBot(5,2)
    endif
    return
end subroutine totfish
