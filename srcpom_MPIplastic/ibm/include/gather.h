
!--maximum number of SIs allowed for every variant (if bigger merge in routine gather)
    integer :: nsiVariants(1:nvariantsmax,nspecies)

    data nsiVariants(1:nvariants1,1)/20000,20000,20000,20000,20000,20000/    ! Micro
    data nsiVariants(1:nvariants2,2)/20000,20000,20000,20000,20000/          ! Macro

    real :: distmin(nvariantsmax),distmax(nvariantsmax)
    data distmin/0.5,0.5,0.5,0.5,0.5,0.5/
    data distmax/200.,200.,200.,200.,200.,200./

