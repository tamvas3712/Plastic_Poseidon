program prepare_netcdf_ibmPlastic_addSources
!!  This program adds plastics sources from Rivers and WWT (micro) and coastal cities/beaches (macro)
    use poseidon_common_subroutines
    use ibm_common_subroutines
    include "netcdf.inc"
    include 'grid.h'
    include 'ibmall.fcm'

    parameter (nsou=3250)!max number of sources

    integer :: nsi_id,species_id,xfi_id,nfi_id,variant_id,xsi_id,ysi_id, &
    ref_id,iarea_id,rbt_id,fou_id,fdep_id,icountdays_id, dista_id

    integer :: h_id,alon_id,alat_id 

    character indf*2
    real :: xxfi,nnfi,xxsi,yysi
    real :: alonsRiv(nsou),alatsRiv(nsou),simicroRiv(nsou),simacroRiv(nsou)
    real :: alonsWWT(nsou),alatsWWT(nsou),simicroWWT(nsou),sinanoWWT(nsou)
    real :: alonsBea(nsou),alatsBea(nsou),simacroBea(nsou)
    integer :: isisRiv(nsou),isis1wwt(nsou),isis2wwt(nsou),isisBea(nsou)
    integer :: itypeWWT(nsou)
    real :: pop(nsou)

    integer :: bot(nmax)
    
    real ::  pload(4),factsize(2,2)
    integer :: NCID1,NCID,status
    character CFILE1*128

!    read(*,110) iday2,imo2 !it reads the date and outputs plastics accordingly.
    write(*,*)iday2,imo2
    110 format(1X,i2,1X,i2)

!---read model grid (used in subroutine init_lagrSI)----------------
   status = NF_OPEN('../in/pom.gridplastic.nc',0, ncid)
    write(*,*)'OK',filename,ncid
!---get variables
    status=nf_inq_varid(ncid,'h',h_id)
    status=nf_inq_varid(ncid,'longitude',alon_id)
    status=nf_inq_varid(ncid,'latitude',alat_id)
    status=nf_get_var_real(ncid,h_id,h)
    status=nf_get_var_real(ncid,alon_id,alon)
    status=nf_get_var_real(ncid,alat_id,alat)
    status = NF_CLOSE (NCID)   
     
    nsi=nmax
    write(*,*)'NSI=',nsi,h(200,50)

!read ibm file
    status = NF_OPEN('../out/pom.ibmplastic.nc',0, ncid)
    write(*,*)'OK',filename,ncid
!---get variables
    status=nf_inq_varid(ncid,'species',species_id)
!    status=nf_inq_varid(ncid,'xfi',xfi_id)
    status=nf_inq_varid(ncid,'nfi',nfi_id)
    status=nf_inq_varid(ncid,'variant',variant_id)
    status=nf_inq_varid(ncid,'xsi',xsi_id)
    status=nf_inq_varid(ncid,'ysi',ysi_id)
    status=nf_inq_varid(ncid,'ref',ref_id)
    status=nf_inq_varid(ncid,'iarea',iarea_id)
    status=nf_inq_varid(ncid,'fdep',fdep_id)
    status=nf_inq_varid(ncid,'rbt',rbt_id)
    status=nf_inq_varid(ncid,'fou',fou_id)
    status=nf_inq_varid(ncid,'icountdays',icountdays_id)
    status=nf_inq_varid(ncid,'dista',dista_id)
    status=nf_get_var_int(ncid,species_id,species)
    status=nf_get_var_int(ncid,variant_id,variant)
    status=nf_get_var_int(ncid,iarea_id,iarea)
    status=nf_get_var_int(ncid,ref_id,ref)
    status=nf_get_var_int(ncid,icountdays_id,icountdays)
!    status=nf_get_var_real(ncid,xfi_id,xfi)
    status=nf_get_var_real(ncid,nfi_id,nfi)
    status=nf_get_var_real(ncid,xsi_id,xsi)
    status=nf_get_var_real(ncid,ysi_id,ysi)
    status=nf_get_var_real(ncid,rbt_id,rbt)
    status=nf_get_var_real(ncid,fou_id,fou)
    status=nf_get_var_real(ncid,fdep_id,fdep)
    status=nf_get_var_real(ncid,dista_id,dista)
    status = NF_CLOSE (NCID)

    !write(*,*)'TEST'

!---find empty SIs to fill with new from sources
    call totfish(1,iempty)

!---random density increase of bottles-----------------------------
    icount=0
    do n=1,nsi
        if(species(n) == 2 .AND. variant(n) == 3 .AND. ref(n) == 0)then
            icount=icount+1
            bot(icount)=n
        endif
    enddo
    do n=1,int(icount*0.025)
!        do n=1,int(icount*0.012)
        nnn=randomf1(1,icount)
        nn=bot(nnn)
        fou(nn)=1300.
        !write(*,*)'BOTTLE',icount,n,nn,nnn,species(nn),variant(nn)
    enddo

!---decrease partics on beach (burial, removal etc)-----------------
    totmacro1=0.
    do k=1,nsi
        nsp=species(k)
        ivariant=variant(k)
        iref=ref(k)
        if(iref == 77 .and. ivariant .ne. 0)then
    !        factBea=burial(ivariant,nsp)*0.75 !decrease in coarse MED10 (increased
    !        beaching)
            !write(*,*)'before Burial',ivariant,nsp
            factBea=burial(ivariant,nsp)
            nfi(k)=nfi(k)-nfi(k)*factBea
            if(ivariant == 1 .AND. nsp == 2 .AND. ref(k) >= 0) &
            totmacro1=totmacro1+nfi(k)
        endif
    enddo
    write(*,*)'TOTMACRO1-BEACHOUT',totmacro1

   
!<<<<<<<<<<<<<<<<<<<<<<<< ADD SOURCES>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!-----RIVER INPUT----------------------------------------------------!!     

    totmacro3=0.
!---read sources
    icount=0
    totload1=0.
    totload2=0.
    write(indf,'(i2.2)')imo2
    open(14,file='../data/sourcesMedEquaModFillAnnNew.dat') 
!! alonsRiv,alatsRiv = source location
!! isisRiv=number of Super Individuals (SI) added
!! simicroRiv,simacroRiv= number of micro/macro particles in each SI
!! korigin = unique number identifying the source
    do k=1,nsou
        read(14,188,end=1000)alonsRiv(k),alatsRiv(k),isisRiv(k),simicroRiv(k),simacroRiv(k),korigin
        icount=icount+1
        188 format(2(f6.2,1x),I4,1x,2(e12.5,1x),I4)
        !write(*,*)alonsRiv(k),alatsRiv(k),isisRiv(k),simicroRiv(k),simacroRiv(k),korigin
        !write(*,*)'SOURCE',k,isisRiv(k)
        xxsi=alonsRiv(k)
        yysi=alatsRiv(k)
        totload1=totload1+isisRiv(k)*simicroRiv(k)*0.00036
        totload2=totload2+isisRiv(k)*simacroRiv(k)*0.17

        do iloop=1,isisRiv(k)   !input of SIs from River

    !---microplastics small (No Input for now) 
    !       nsize=2
    !       do isize=1,nsize
    !       nsp=1
    !       ivariant=isize
    !       xxfi=4/3*3.14*(0.5*diam(ivariant,1))**3.*dens(ivariant,1)*1000 !weight(gr)(m=V*rho, V=4/3*R^3)
    !       nnfi=simicro(k)*0.7*10.
    !       nnfi=simicro(k)*0.
    !       call StartNew(nsp,ivariant,xxfi,nnfi,xxsi,yysi,korigin) !micro
    !   enddo

        !---microplastics large
            scalefactor=1.3
            ttload=0
            nsize=4
            do isize=1,nsize
                ivariant=2+isize
                pload(isize)=1./(diam(ivariant,1)*0.5)**scalefactor
                ttload=ttload+pload(isize)
            enddo

            do isize=1,nsize
                pload(isize)=pload(isize)/ttload
            enddo

            do isize=1,nsize
                nsp=1
                ivariant=2+isize
                xxfi=4/3*3.14*(0.5*diam(ivariant,1))**3.* &
                dens(ivariant,1)*1000 !weight(gr)(m=V*rho,V=4/3*R^3)
                nnfi=simicroRiv(k)*pload(isize)
                call StartNew(nsp,ivariant,xxfi,nnfi,xxsi,yysi,korigin) !nano
            enddo

        !--macroplastics
           factSou=1.1
           nsp=2
            do ivariant=1,5
               call random_number(r)
               xxfi=xfiInit(ivariant,nsp)
               nnfi=simacroRiv(k)*riverMacro(ivariant)*riverMacro20(ivariant)*factSou
               call StartNew(nsp,ivariant,xxfi,nnfi,xxsi,yysi,korigin)
            enddo

        enddo
    enddo
    1000 continue
    close(14)

    call totfish(1,iempty)
    write(*,*)'River Sources added'
    write(*,*)'TOTAL BIOMASS MICRO/MACRO',365*totload1*1.e-6, &
    365*totload2*1.e-6

    write(*,*)'TOTMACRO1-RIVIN',totmacro1+totmacro2,totmacro3, &
    totmacro1+totmacro2+totmacro3

    if(iempty < 50000)then
        call gathersiPlastic
    endif

!-----WWT INPUT (ONLY MICRO)--------------------------------------------!! 

!---read sources
    icount=0
    totsis=0
    totload1=0.
    totload2=0.
    totpart1=0.
    totpart2=0.

    open(124,file='../data/wwtMedpoptype.dat')
    do i=1,nsou
        read(124,*,end=1005)pop(i),icountry,itypeWWT(i),ratiotreat
        if(itypeWWT(i) < 0.)itypeWWT(i)=1
    enddo
    close(124)
  1005 continue

    open(14,file='../data/sourcesMed_wwt.dat')
!    open(15,file='../data/sourcesMedWWT.csv')
!! alonsWWT,alatsWWT = source location
!! isis1WWT=number of Super Individuals (SI) added for small micro (<300um)
!! isis2WWT=number of Super Individuals (SI) added for large micro (>300um)
!! simicroWWT,sinanoWWT= number of large/small micro particles in each SI
!! korigin = unique number identifying the source
    do k=1,nsou                 !read sources
        read(14,189,end=1010)alonsWWT(k),alatsWWT(k), &
        isis2WWT(k),isis1WWT(k),simicroWWT(k),sinanoWWT(k),korigin
        icount=icount+1

        write(*,*)'SOURCE',k,isis1WWT(k),isis2WWT(k)
        totsis=totsis+isis1WWT(k)+isis2WWT(k)
        189 format(2(f6.2,1x),I4,1x,I4,1x,2(e12.5,1x),I4)

        xxsi=alonsWWT(k)
        yysi=alatsWWT(k)

!  input of small micro SIs from WWT
        nsize=2
        factsize(1,1)=0.7
        factsize(2,1)=0.3
        factsize(1,2)=0.5
        factsize(2,2)=0.5
        fact=1

        do iloop=1,isis1WWT(k) 
            do isize=1,nsize
                if(itypeWWT(korigin-5000) > 1)then !assume 70% = (20–100 µm) when treated (Talvitie et al., 2017), otherwise 50%
                    fact=factsize(isize,1)     
                else
                    fact=factsize(isize,2)
                endif
                nsp=1
                ivariant=isize
                xxfi=4/3*3.14*(0.5*diam(ivariant,1))**3.* &
                dens(ivariant,1)*1000       !weight(gr)(m=V*rho, V=4/3*R^3)
                nnfi=sinanoWWT(k)*fact
                call StartNew(nsp,ivariant,xxfi,nnfi,xxsi,yysi,korigin) !micro
                totload1=totload1+nnfi*xxfi
                totpart1=totpart1+nnfi
            enddo
        enddo

!  input of large (>300um) micro SIs from WWT (only untreated waste) / assume  #particles=f(1/size)^scalefactor
        nsize=4
        scalefactor=1.3
        ttload=0
        do isize=1,nsize
            ivariant=2+isize
            pload(isize)=1./(diam(ivariant,1)*0.5)**scalefactor
            ttload=ttload+pload(isize)
        enddo
        do isize=1,nsize
            pload(isize)=pload(isize)/ttload !normalize with total load
        enddo

        do iloop=1,isis2WWT(k)
            do isize=1,nsize
                nsp=1
                ivariant=2+isize
                xxfi=4/3*3.14*(0.5*diam(ivariant,1))**3.* &
                dens(ivariant,1)*1000 !weight(gr)(m=V*rho, V=4/3*R^3)

                nnfi=simicroWWT(k)*pload(isize)
                call StartNew(nsp,ivariant,xxfi,nnfi,xxsi,yysi,korigin) !nano
                totload2=totload2+nnfi*xxfi
                totpart2=totpart2+nnfi
            enddo
        enddo
        678 continue
    enddo

    1010 continue

    close(14)
    call totfish(1,iempty)

    write(*,*)'WWT Sources added',totsis,icount
    write(*,*)'TOTAL PARTICLES<300/>300um',totpart1,totpart2
    write(*,*)'TOTAL BIOMASS<300/>300um',365*totload1*1.e-6, &
    365*totload2*1.e-6


!---COASTAL CITIES/BEACHES INPUT (ONLY MACRO)---------------------------------!!
    totmacro=0.
    icount=0
    totsis=0.
    open(14,file='../data/sourcesMed_coast.dat')
!! alonsBea,alatsBea = source location
!! isisBea=number of Super Individuals (SI) added
!! simacroBea= number of macro particles in each SI
!! korigin = unique number identifying the source
    do k=1,nsou
        read(14,190,end=1020)alonsBea(k),alatsBea(k), &  !---read sources
        isisBea(k),simacroBea(k),korigin
        icount=icount+1
        write(*,*)'SOURCE',k,isisBea(k)
        190 format(2(f6.2,1x),I4,1x,e12.5,1x,I4)

        xxsi=alonsBea(k)
        yysi=alatsBea(k)
        factSou=1.

        totsis=totsis+isisBea(k)

       do iloop=1,isisBea(k)
         nsp=2
        do ivariant=1,5
         xxfi=xfiInit(ivariant,nsp)
         nnfi=simacroBea(k)*beachMacro(ivariant)*beachMacro20(ivariant)*factSou
         call StartNew(nsp,ivariant,xxfi,nnfi,xxsi,yysi,korigin)
         totmacro=totmacro+nnfi*xxfi
        enddo
       enddo

    enddo
    1020 continue
    close(14)
    call totfish(1,iempty)

    write(*,*)'Beach Sources added',totsis,icount
    write(*,*)'TOTAL BIOMASS BEACH MACRO',365*totmacro*1.e-6

    call gathersiPlastic

    call totfish(1,iempty)

    if(iempty < 50000)then !repeat gather if empty<50000
        call gathersiPlastic
    endif

!---open NetCDF files
    CFILE1='../data/pom.ibmplasticNew.nc'
    status=NF_CREATE(CFILE1(1:(INDEX(CFILE1,' '))-1),NF_SHARE,NCID1)!create dataset
    write(*,*)status,NCID1
    status = NF_REDEF(NCID1)
    status = NF_DEF_DIM(NCID1, 'nsi',nsi, nsi_id)
    write(*,*)status,nsi,nsi_id,nsi
    status = NF_REDEF(NCID1)
    status = NF_DEF_VAR(NCID1,'species',NF_INT,1,nsi_id,species_id)
!    status = NF_DEF_VAR(NCID1,'xfi',NF_FLOAT,1,nsi_id,xfi_id)
    status = NF_DEF_VAR(NCID1,'nfi',NF_FLOAT,1,nsi_id,nfi_id)
    status = NF_DEF_VAR(NCID1,'variant',NF_INT,1,nsi_id,variant_id)
    status = NF_DEF_VAR(NCID1,'xsi',NF_FLOAT,1,nsi_id,xsi_id)
    status = NF_DEF_VAR(NCID1,'ysi',NF_FLOAT,1,nsi_id,ysi_id)
    status = NF_DEF_VAR(NCID1,'ref',NF_INT,1,nsi_id,ref_id)
    status = NF_DEF_VAR(NCID1,'iarea',NF_INT,1,nsi_id,iarea_id)
    status = NF_DEF_VAR(NCID1,'fdep',NF_FLOAT,1,nsi_id,fdep_id)
    status = NF_DEF_VAR(NCID1,'rbt',NF_FLOAT,1,nsi_id,rbt_id)
    status = NF_DEF_VAR(NCID1,'fou',NF_FLOAT,1,nsi_id,fou_id)
    status = NF_DEF_VAR(NCID1,'dista',NF_FLOAT,1,nsi_id,dista_id)
    status = NF_DEF_VAR(NCID1,'icountdays',NF_INT,1,nsi_id, &
    icountdays_id)
    status = NF_ENDDEF(ncid1, iret)
    write(*,*)CFILE1
    status = nf_put_var_int(ncid1,species_id,species)
!    status = nf_put_var_real(ncid1,xfi_id,xfi)
    status = nf_put_var_real(ncid1,nfi_id,nfi)
    status = nf_put_var_int(ncid1,variant_id,variant)
    status = nf_put_var_real(ncid1,xsi_id,xsi)
    status = nf_put_var_real(ncid1,ysi_id,ysi)
    status = nf_put_var_int(ncid1,ref_id,ref)
    status = nf_put_var_int(ncid1,iarea_id,iarea)
    status = nf_put_var_int(ncid1,icountdays_id,icountdays)
    status = nf_put_var_real(ncid1,fdep_id,fdep)
    status = nf_put_var_real(ncid1,rbt_id,rbt)
    status = nf_put_var_real(ncid1,fou_id,fou)
    status = nf_put_var_real(ncid1,dista_id,dista)
    status = NF_CLOSE (NCID1)
    stop
end program







subroutine StartNew(nsp,ivariant,xxfi,nnfi,xxsi,yysi,korigin)
!! ??
    use ibm_common_subroutines
    implicit none
    include 'ibmall.fcm'
    integer :: nnn,nn,n,iflagFill,nsp,ivariant,korigin
    real :: xxfi,nnfi,xxsi,yysi,xxsiNew,yysiNew

!---first fill empty S.I
    iflagFill=0
    do nn=1,nmax
        if(nempty(nn) /= 0)then
            nnn=nempty(nn)
            species(nnn)=nsp
            variant(nnn)=ivariant
            iflagFill=1
            nfi(nnn)=nnfi
            xfi(nnn)=xxfi
            ref(nnn)=0
            iarea(nnn)=korigin
            icountdays(nnn)=0
            icountdays(nnn)=0
            dista(nnn)=0.
            if(nsp == 2 .AND. ivariant <= 4)then
                fou(nnn)=dens(ivariant,nsp)
            else
                fou(nnn)=0.
            endif
            fdep(nnn)=0.
            rbt(nnn)=0.
            nempty(nn)=0
!            write(*,*)'NEW SI ON EMPTY',nnn,nsp,ivariant,nnfi
            goto 145
        endif
    enddo

    145 continue
    if(iflagFill == 0)then
        write(*,*)'ERROR EMPTY'
    else
        call init_lagrSI(xxsi,yysi,xxsiNew, yysiNew)

        xsi(nnn)=xxsiNew
        ysi(nnn)=yysiNew
!        write(*,*)'OK NEW SI',xxsi,yysi,xxsiNew, yysiNew,ivariant
!        write(*,*)'OK NEW SI',nsp,ivariant,nnfi,nnn
    endif
    return
end subroutine StartNew






subroutine totfish(iwrite,icountempty)
!! ??
    include 'ibmall.fcm'
    integer :: ch,i,nn,iwres
    integer :: icountempty,iwrite
    real :: weira(nvariantsmax,nspecies)
    real :: countd(nvariantsmax,nspecies)

    do nsp=1,nspecies
        do i=1,nvariants(nsp)
            weivariant(i,nsp)=0.
            tbvariant(i,nsp)=0.
            tnvariant(i,nsp)=0.
            tsivariant(i,nsp)=0.
            countd(i,nsp)=0.
            tsivariantBea(i,nsp)=0.
            tsivariantBot(i,nsp)=0.
        enddo
    enddo

!---Calculation of empty SIs,then fill with new eggs in ersem_fishGrad.F
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
            tbvariant(ivariant,nsp)=tbvariant(ivariant,nsp)+nfi(n)*xfiInit(ivariant,nsp)
            tnvariant(ivariant,nsp)=tnvariant(ivariant,nsp)+nfi(n)
            tsivariant(ivariant,nsp)=tsivariant(ivariant,nsp)+1
            if(ref(n) > 0)tsivariantBea(i,nsp)=tsivariantBea(i,nsp)+1
            if(ref(n) == -1)tsivariantBot(i,nsp)=tsivariantBot(i,nsp)+1
            weivariant(ivariant,nsp)=weivariant(ivariant,nsp)+nfi(n)*xfiInit(ivariant,nsp)
            countd(ivariant,nsp)=countd(ivariant,nsp)+icountdays(n)*nfi(n)
        endif
    enddo

    do n=1,nspecies
        do i=1,nvariants(n)
            if(tnvariant(i,n) > 0.)then
                weivariant(i,n)=weivariant(i,n)/tnvariant(i,n)
                countd(i,n)=countd(i,n)/tnvariant(i,n)
            else
                weivariant(i,n)=0.
                countd(i,n)=0.
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
                write(*,*)'TOTAL BIOMASS',ivariant,tbvariant(ivariant,n)*1.e-6
                write(*,*)'PARTICLES',ivariant,tnvariant(ivariant,n)
                write(*,*)'SIs',ivariant,tsivariant(ivariant,n)
                write(*,*)'SIs-BEACH',ivariant,tsivariantBea(ivariant,n)
                write(*,*)'SIs-BOT',ivariant,tsivariantBot(ivariant,n)
!                write(*,*)'PARTICLE WEIGHT(mgr)',ivariant,weivariant(ivariant,n)
                write(*,*)'DAYS',ivariant,countd(ivariant,n)
            enddo
        enddo
    endif
    return
end subroutine totfish







subroutine init_lagrSI(xilonSI,yilatSI,xxsi, yysi)
!! initialize lagrangian Super Individuals.
    use ibm_common_subroutines
    include 'ibmall.fcm'   
    include 'grid.h'
    include 'lagr.h'

    real :: xi,yi !position in cartesian coordinates (m)
    real :: xn1,yn1,zn1
    real ::     pi,                g
    parameter ( pi = 3.1415926536, g = 9.81 )


!---Proceed with Initialization - Get Geographic coordinates
!    and Set Cartesian and Grid coordinates
    call Geo2Cart(xi, yi, alon(1,1)-(stepx/2.), &
        alat(1,1)-(stepy/2.), xilonSI,yilatSI, pi, dg)

!     write(*,*)'OK1-initLagr',xilonSI,yilatSI,xxsi, yysi,xi,yi

!---position of each particle
    50 continue
    r1 = rand(0)
    r2 = rand(0)
    xn1 = xi + r1 * rin * cos(2. * pi * r2)
    yn1 = yi + r1 * rin * sin(2. * pi * r2)
!    zn = 10.
    zn1 = 0.

    call FCoord(io, jo, alon(1,1)-(stepx/2.), &
    alat(1,1)-(stepy/2.), stepx, stepy, &
    xn1, yn1, pi, dg)

!    write(*,*)'TEST',io, jo,yn,yi,r2,sin(2. * pi * r2),pi

    if(io < 1 .OR. io > im .OR. jo < 1 .OR. jo > jm) goto 50
    if(h(io, jo) <= 1.) goto 50

!---initialize ibm (i,j) and position
    call Cart2Geo(xxsi, yysi, alon(1,1)-(stepx/2.), &
    alat(1,1)-(stepy/2.), xn1, yn1, pi, dg)

    if(abs(xxsi-xilonSI) > 1.) &
    write(*,*)'RANDX',xi,xn1,xi-xn1,r1,r2,rin
    if(abs(yysi-yilatSI) > 1.) &
    write(*,*)'RANDY',yi,yn1,yi-yn1,r1,r2,rin
    return
end subroutine init_lagrSI



