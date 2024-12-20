program write_netcdf_ibmPlastic
!! This program reads INIT.DAT and writes netcdf ibm plastic data
    use poseidon_common_subroutines,only:findtime0
    include "netcdf.inc"
    include 'grid.h'
    include 'ibmall.fcm'

    parameter (nsi_local = 3000,nsou=3250)
    integer :: ch
    integer :: nsi_id,species_id,xfi_id,nfi_id,variant_id,xsi_id,ysi_id, &
    iarea_id,ref_id,fdep_id,icountdays_id,fou_id,dista_id
    integer :: NCID1,NCID,status
    CHARACTER CFILE1*128,indf1*2,indf2*2,indf3*2,indf4*2,filename*20
     
!--read INIT.DAT (grid+TSclim,Rmean,TS at boundaries)
    OPEN(50,FILE='../data/INITnew.DAT',FORM='UNFORMATTED')
    READ(50) Z,ZZ,DZ,DZZ,ALON,ALAT,DX,DY,H
    CLOSE(50)

    nsi=nmax
    write(*,*)'NSI=',nsi
      
    status = NF_OPEN('../in/pom.ibmplastic.nc',0, ncid)

    write(*,*)'OK',filename,ncid
! get variables
    status=nf_inq_varid(ncid,'species',species_id)
    status=nf_inq_varid(ncid,'xfi',xfi_id)
    status=nf_inq_varid(ncid,'nfi',nfi_id)
    status=nf_inq_varid(ncid,'variant',variant_id)
    status=nf_inq_varid(ncid,'xsi',xsi_id)
    status=nf_inq_varid(ncid,'ysi',ysi_id)
    status=nf_inq_varid(ncid,'ref',ref_id)
    status=nf_inq_varid(ncid,'iarea',iarea_id)
    status=nf_inq_varid(ncid,'icountdays',icountdays_id)
    status=nf_inq_varid(ncid,'fdep',fdep_id)
    status=nf_inq_varid(ncid,'fou',fou_id)
    status=nf_inq_varid(ncid,'dista',dista_id)
    status=nf_get_var_int(ncid,species_id,species)
    status=nf_get_var_int(ncid,variant_id,variant)
    status=nf_get_var_int(ncid,iarea_id,iarea)
    status=nf_get_var_int(ncid,ref_id,ref)
    status=nf_get_var_int(ncid,icountdays_id,icountdays)
    status=nf_get_var_real(ncid,xfi_id,xfi)
    status=nf_get_var_real(ncid,nfi_id,nfi)
    status=nf_get_var_real(ncid,xsi_id,xsi)
    status=nf_get_var_real(ncid,ysi_id,ysi)
    status=nf_get_var_real(ncid,fdep_id,fdep)
    status=nf_get_var_real(ncid,fou_id,fou)
    status=nf_get_var_real(ncid,dista_id,dista)
    status = NF_CLOSE (NCID)

    call totfish(1)

!--write binary file
    read(*,110) iday2,imo2,iyr2
    write(*,*)iday2,imo2,iyr2
    110 format(1X,i2,1X,i2,1X,i2)

    iyr=iyr2
    if(iyr2 <= 50)then
        iyr=iyr+2000
    else
        iyr=iyr+1900
    endif

    write(indf1,'(i2.2)')iday2
    write(indf2,'(i2.2)')imo2
    write(indf3,'(i2.2)')iyr2

    open(99,file='FISH'//indf3//indf2//indf1//'.18UTC',form='unformatted')
    !write(99)nsi
    do n=1,nsi
        if(variant(n) > 0)then
            write(99)n
            write(99)species(n),xfi(n),nfi(n),variant(n),xsi(n),ysi(n),ref(n),fdep(n), &
            icountdays(n),fou(n),dista(n)
        endif
    enddo
    close(99)

end program      







subroutine totfish(iwrite)
!! ??
    include 'ibmall.fcm'
    integer :: ch,i,nn,iwres
    integer :: icountempty,iwrite
    real :: weira(nvariantsmax,nspecies)

    do nsp=1,nspecies
        do i=1,nvariants(nsp)
            weivariant(i,nsp)=0.
            tbvariant(i,nsp)=0.
            tnvariant(i,nsp)=0.
            tsivariant(i,nsp)=0.
        enddo
    enddo

! Calculation of empty SIs,then fill with new eggs in ersem_fishGrad.F
    do n=1,nmax
        nempty(n)=0
    enddo
    icountempty=0
    icountsource=0
    icountback=0
    icountbea=0
    icountbot=0
    totbiom1=0.
    totbiom2=0.
    totbiom1bea=0.
    totbiom2bea=0.

    do n=1,nsi
        if(variant(n) == 0)then
            icountempty=icountempty+1
            nempty(icountempty)=n
        else
            if(iarea(n) > 1)icountsource=icountsource+1
            if(iarea(n) == 1)icountback=icountback+1
            if(ref(n) == 77)icountbea=icountbea+1
            if(ref(n) == -1)icountbot=icountbot+1
            nsp=species(n)
            ivariant=variant(n)
            tbvariant(ivariant,nsp)=tbvariant(ivariant,nsp)+ &
            nfi(n)*xfiInit(ivariant,nsp)
            if(nsp == 1)totbiom1=totbiom1+nfi(n)*xfiInit(ivariant,nsp)
            if(nsp == 2)totbiom2=totbiom2+nfi(n)*xfiInit(ivariant,nsp)
            if(nsp == 1 .AND. ref(n) == 77)totbiom1bea=totbiom1bea+ &
            nfi(n)*xfiInit(ivariant,nsp)
            if(nsp == 2 .AND. ref(n) == 77)totbiom2bea=totbiom2bea+ &
            nfi(n)*xfiInit(ivariant,nsp)
            tnvariant(ivariant,nsp)=tnvariant(ivariant,nsp)+nfi(n)
            tsivariant(ivariant,nsp)=tsivariant(ivariant,nsp)+1
            weivariant(ivariant,nsp)=weivariant(ivariant,nsp)+ &
            nfi(n)*xfiInit(ivariant,nsp)
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
        write(*,*)'SIs Source',icountsource
        write(*,*)'SIs Background',icountback
        write(*,*)'SIs Beach',icountbea
        write(*,*)'SIs BOTTOM',icountbot
        do n=1,nspecies
            write(*,*)'SPECIES',n
            adultbiom=0.
            do ivariant=1,nvariants(n)
                write(*,*)'BIOMASS',ivariant,tbvariant(ivariant,n)*1.e-6
                write(*,*)'PARTICLES',ivariant,tnvariant(ivariant,n)
                write(*,*)'SIs',ivariant,tsivariant(ivariant,n)
!                write(*,*)'PARTICLE WEIGHT(mgr)',ivariant,weivariant(ivariant,n)
            enddo
        enddo
    endif
    write(*,*)'BIOMASS MICROSEA/BEACH',totbiom1*1.e-6, &
    totbiom1bea*1.e-6
    write(*,*)'BIOMASS MACRO SEA/BEACH',totbiom2*1.e-6, &
    totbiom2bea*1.e-6
    return
end subroutine totfish





