subroutine gathersiPlastic
!! In order to avoid a high increase of Super Individuals (SIs), they are merged together if they are proximal
    use poseidon_common_subroutines
    use ibm_common_subroutines
    include 'ibmall.fcm'
    include 'grid.h'
    include 'lagr.h'
    include 'gather.h'

    integer :: maxsi
    parameter (maxsi=10000,igrid=4)
    parameter (im2=im/igrid+2,jm2=jm/igrid+2) !was +1 before, compiler flooring?
    integer ::  siloc(im2,jm2,maxsi)
    integer ::  nsiloc(im2,jm2)
    real ::     ansiloc(im2,jm2)
    integer ::  idone(im2*jm2)
    integer :: variant1,variant2,icountday1,icountday2
    real :: nfi1,nfi2,dista1,dista2
    integer :: nsivariant(nvariantsmax,nspecies),nsivariantall(nvariantsmax,nspecies)
    integer :: ngath(nvariantsmax,nspecies),gathall(nvariantsmax,nspecies)
    integer :: iir(3)
    integer :: iivariant,k,n,ivariant,n1,n2,nsp,iref,iiref,ref1,ref2
    real :: fact(3)
    data iir/0,77,-1/
    data fact/0.6,0.2,0.1/
    write(*,*)'START GATHER'

    xs=(alon(2,1)-alon(1,1))*igrid
    ys=(alat(1,2)-alat(1,1))*igrid

    do 5 nnsp=1,nspecies
 !------only microplastics (merge background with sources)
        do 10 iivariant=1,nvariants(nnsp)
            gathall(iivariant,nnsp)=0
            nsivariantall(iivariant,nnsp)=0
            do 12 ir=1,3
                iiref=iir(ir)
                ffact=fact(ir)
                write(*,*)'GATH',nnsp,iivariant,iiref
                nsivariant(iivariant,nnsp)=0
                ngath(iivariant,nnsp)=0
                do k=1,maxsi
                    do j=1,jm2
                        do i=1,im2
                            siloc(i,j,k)=0
                        enddo
                    enddo
                enddo
                do j=1,jm2
                    do i=1,im2
                        nsiloc(i,j)=0
                    enddo
                enddo

!-------------- place all si in the 2-D grid
                do n=1,nsi
                    x=xsi(n)
                    y=ysi(n)
                    ii=int((x-(alon(1,1)-xs*0.5))/xs)+1
                    jj=int((y-(alat(1,1)-ys*0.5))/ys)+1
                    ivariant=variant(n)
                    nsp=species(n)
                    iiarea=iarea(n)
                    iref=ref(n)
                    if(ivariant .eq. iivariant .and. nsp .eq. nnsp &
                !     .and..iiarea.gt.1.and.iref.eq.iiref)then !gather only SIs from sources (keep initial SIs)
                    .and. iref .eq. iiref)then
                !        if(variant(n).eq.iivariant.and.nsp.eq.nnsp)then
                        nsivariant(ivariant,nnsp)=nsivariant(ivariant,nnsp)+1 !calculate #of SIs for each variant
                        nsivariantall(ivariant,nnsp)=nsivariantall(ivariant,nnsp)+1

                        if(variant(n)==0)write(*,*)'IVARIANT-ZERO',n

                        do k=1,maxsi !why ??? tamvas
                            if(siloc(ii,jj,k) .eq. 0)then
                                siloc(ii,jj,k)=n
                                nsiloc(ii,jj)=k
                                ansiloc(ii,jj)=float(nsiloc(ii,jj))
                                goto 123
                            endif
                        enddo
                        123 continue

                !        if(nsiloc(ii,jj).lt.maxsi.and.siloc(ii,jj,k).gt.0)then
                    !        nsiloc(ii,jj)=nsiloc(ii,jj)+1
                    !        k=nsiloc(ii,jj)
                    !        siloc(ii,jj,k)=n
                !        endif

                    endif
                enddo

                write(*,*)'OK',nnsp,iivariant,nsivariant(iivariant,nnsp),nsiVariants(iivariant,nnsp)

                if(nsivariant(iivariant,nnsp) .le. nsiVariants(iivariant,nnsp)*ffact)goto 12

        !        write(*,*)'NSI',nsi
        !        CALL PRXY('NSILOC',TIME,ansiloc(1,1),IM,1,JM,1,1.E+1)
        !        stop

 !------------- gather si when in the same grid point within specified radius
                slon=alon(1,1)-(xs/2.)
                slat=alat(1,1)-(ys/2.)
                ngathall=0
                do n=1,im2*jm2
                    j=(n-1)/im2+1
                    i=mod(n-1,im2)+1
                    if(nsiloc(i,j) .gt. 1)then
                        ngathall=ngathall+1
                        idone(ngathall)=n
                    endif
                enddo
                write(*,*)'GATHER',iiref,ngathall,nsiVariants(iivariant,nnsp)*ffact,nsivariant(iivariant,nnsp)

              !  write(*,*)'TEST-SILOC',siloc(42,30,7),variant(siloc(42,30,7))

                do 15 iloop=1,500
                    do 20 nn=1,ngathall
                        nnn=randomf(1,ngathall)
                        n=idone(nnn)
                        j=(n-1)/im2+1
                        i=mod(n-1,im2)+1
                !        do k=1,nsiloc(i,j)-1
                        do 25 k=1,1
                            n1=siloc(i,j,k)
                            species1=species(n1)
                            variant1=variant(n1)
                            !if(variant1==0)write(*,*)'IVARIANT1-ZERO',n1,i,j,k
                            icountday1=icountdays(n1)
                            nfi1=nfi(n1)
                            dista1=dista(n1)
                            iarea1=iarea(n1)
                            ref1=ref(n1)
                    !        if(nnsp.eq.1.and.iivariant.eq.1)write(*,*)'TEST1',ref1,nn,k
                         !   write(*,*)'TEST',i,j,k,nsiloc(i,j)
                     !----- calculate x,y coordinates in (m)
                            call Geo2Cart(x1, y1, slon, slat, xsi(n1),ysi(n1), pi, dg)
                         !   write(*,*)'TEST2',i,j,k
                            do 30 kk=k+1,nsiloc(i,j)
                                n2=siloc(i,j,kk)
                                species2=species(n2)
                                variant2=variant(n2)
                                if(variant2.eq.0)goto 30
                          !      if(variant2.eq.0)write(*,*)'ZERO,should exit',n2,i,j,kk
                                icountday2=icountdays(n2)
                                nfi2=nfi(n2)
                                dista2=dista(n2)
                                iarea2=iarea(n2)
                                ref2=ref(n2)
                                call Geo2Cart(x2, y2, slon, slat, xsi(n2),ysi(n2), pi, dg)
                                dist=sqrt((x1-x2)**2.+(y1-y2)**2.)*1.E-3 ! distance in Km
                          !      write(*,*)'TEST3',i,j,k,kk,nn,variant2,n2
                     !           if(variant2.eq.0)then !test
                     !               ddist=1!test
                     !           else !test
                           !         ddist=distmax(variant2) !problematic if variant2=0
                                if(variant2.ne.0)ddist=distmax(variant2)
                     !           endif !test
                                !if(variant2.eq.0)write(*,*)'ZEROAGE',kk,n2
                                
                    !            if(nnsp.eq.1.and.iivariant.eq.1.and.ir.eq.2)write(*,*)'TEST2',i,j,kk,ref1,ref2,dist
                                if(species1 .eq. species2 .and. variant1 .eq. variant2 .and. dist .lt. ddist &
                                 .and. iarea2*iarea1 .gt. 1 &
                                 .and. ref1 .eq. ref2) then
                            !         if(iarea1+iarea2.gt.2)then  !gather sources with background (not background with background)
                                    if(variant2.eq.0.or.variant1.eq.0)write(*,*)'ZEROVARIANT',n1,n2 

                                    if(nsivariant(variant1,nnsp) .gt. nsiVariants(variant1,nnsp)*ffact)then
                                        nsivariant(variant1,nnsp)=nsivariant(variant1,nnsp)-1
                                        nsivariantall(variant1,nnsp)=nsivariantall(variant1,nnsp)-1
                                        ngath(variant1,nnsp)=ngath(variant1,nnsp)+1
                                !        write(*,*)'GATHER',i,j,k,kk,variant1,nsivariant(variant1,nnsp),ngath(variant1,nnsp)
                                        if(iarea1 .eq. 1 .and. iarea2 .gt. 1)iarea(n1)=iarea2
                                        variant(n2)=0
                                        gathered(n2)=n1
                                        if(ref1.eq.0 .or. ref1 .eq. -1)then
                                            xsi(n1)=0.5*(xsi(n1)+xsi(n2))
                                            ysi(n1)=0.5*(ysi(n1)+ysi(n2))
                                        endif
                                        xfi(n1)=(xfi(n1)*real(nfi(n1))+xfi(n2)*real(nfi(n2)))/(nfi(n1)+nfi(n2))
                                        fdep(n1)=(fdep(n1)*real(nfi(n1))+fdep(n2)*real(nfi(n2)))/(nfi(n1)+nfi(n2))
                                        fou(n1)=(fou(n1)*real(nfi(n1))+fou(n2)*real(nfi(n2)))/(nfi(n1)+nfi(n2))
                                        dista(n1)=(dista(n1)*real(dista(n1))+dista(n2)*real(nfi(n2)))/(nfi(n1)+nfi(n2))
                                        if(dista(n1) .gt. 1.e+35)dista(n1)=0.
                                        nfi(n1)=nfi(n1)+nfi(n2)
                                        icountdays(n1)=(icountdays(n1)*real(nfi(n1))+icountdays(n2)*real(nfi(n2)))/ &
                                        (nfi(n1)+nfi(n2))
                                    else
                                        goto 125
                                    endif
                                endif
                            30 continue
                        25 continue
                    20 continue
                15 continue
                125 continue

                if(ngath(iivariant,nnsp) .gt. 0)then
                    write(*,*)'GATHERED',ir,nnsp,iivariant,ngath(iivariant,nnsp),nsivariant(iivariant,nnsp)
                    gathall(iivariant,nnsp)=gathall(iivariant,nnsp)+ngath(iivariant,nnsp)
                endif
            12 continue
        10 continue
!        CALL PRXY('NSILOC',TIME,ansiloc(1,1),IM,1,JM,1,1.E+1)
    5 continue

    do nsp=1,nspecies
        do i=1,nvariants(nsp)
            if(gathall(i,nsp) > 0)then
                write(*,*)'GATHERED-ALL',nsp,i,gathall(i,nsp),nsivariantall(i,nsp)
            endif
        enddo
    enddo
    return
end subroutine gathersiPlastic

