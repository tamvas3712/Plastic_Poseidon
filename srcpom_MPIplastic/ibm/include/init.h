!--initial number of SIs per cell
    integer :: nsiInit(nvariantsmax,nspecies)        !Initial number of Super Individuals
    data nsiInit(1:nvariants1,1)/1,1,1,1,1,1/        !Micro
    data nsiInit(1:nvariants2,2)/1,1,1,1,1/          !Macro

    integer :: icvariant
    common /blkibm/ icvariant(nvariantsmax,nspecies)     !used for initial SI numbering


