! --- Univ. of Athens - Ocean
program pre_meteo_posII
!! Prepare atmospheric data to be fed to the main program. 
    include 'meteo_posII.h'
    parameter (NX=NXP , NY = NYP)
    parameter (NFIELDS = 9)
    character cdate*17, cvar*5
    dimension aux(nx,ny),store_1(nx,ny,nfields)
!---fields file names
    character indf*3,fname1*11,chara*1,directory_fields*80
    character fname*100
    character indf1*2,indf2*2,indf3*2

    integer :: HOURS_TO_INTEGRATE
    data HOURS_TO_INTEGRATE/24/

    DIRECTORY_FIELDS='/data/meteo2/POS_2/'
    DATA FNAME1/'.18/FIELDS.'/

    ldf = 0
    do ii = 1,80
        chara = directory_fields(ii:ii)
        if (chara == '') goto 10
        ldf = ldf + 1
    end do
    10 continue

!---open output files: FIELDS.INPUT and RADIATION.INPUT
    OPEN(50,file='../data/FIELDS.INPUT',form='unformatted')

    read(*,110) iday2,imo2,iyr2
    write(*,*)iday2,imo2,iyr2
    110 format(1X,i2,1X,i2,1X,i2)

    write(indf1,'(i2.2)')iday2
    write(indf2,'(i2.2)')imo2
    write(indf3,'(i2.2)')iyr2

! ----------------------------------------------------
    last_f = hours_to_integrate + 1

!---do first the fields files (tair_2m,u_10,v_10,rh_2m,rain)
    index_f = 6
    do 500 nsteps = 1, last_f,3
        write(indf,'(i3)') index_f
        if (index_f < 10) indf(1:2) = '00'
        if (index_f < 100) indf(1:1) = '0'
        fname =directory_fields(1:ldf)//indf1//indf2//indf3// &
        fname1//indf
        write(*,*) 'OPENING FILE:',fname
        open(40,file=fname,form='unformatted', &
        status='old',err=1000)
        do 450 k = 1,nfields
            call read2d_40(cdate,nx,ny,aux,cvar,iproc)
            write(*,*)'read',k
        !---if reading for any reason is not correct ...exit
            if (iproc == 0) goto 1000

            do i = 1,nx
                do j = 1,ny
                    store_1(i,j,k) = aux(i,j)
                enddo
            enddo
        450 enddo

        call write2d_50(nx,ny,nfields,store_1)

        if ( nsteps == last_f) then
        !---repeat last fields
            call write2d_50(nx,ny,nfields,store_1)
        end if
        close(40)
        index_f = index_f + 3
    500 enddo
    close(50)
    1000 continue
    stop
end program






subroutine read2d_40(cdate,nx,ny,arr,cvar,iproc)
    dimension ARR(NX,NY)
    character CDATE*17, CVAR*5
    iproc = 0
    read(40,err=800) cdate,cvar
    read(40,err=800) ((arr(i,j), i = 1, nx),j = 1, ny)
    iproc = 1
    800 continue
    return
end subroutine read2d_40





subroutine write2d_50(nx,ny,nf,arr)
    dimension ARR(NX,NY,NF)
    write(50) ((arr(i,j,1), i = 1, nx),j = 1, ny)
    write(50) ((arr(i,j,2), i = 1, nx),j = 1, ny)
    write(50) ((arr(i,j,3), i = 1, nx),j = 1, ny)
    write(50) ((arr(i,j,4), i = 1, nx),j = 1, ny)
    write(50) ((arr(i,j,5), i = 1, nx),j = 1, ny)
    write(50) ((arr(i,j,7), i = 1, nx),j = 1, ny)
    write(50) ((arr(i,j,8), i = 1, nx),j = 1, ny)
    return
end subroutine write2d_50
