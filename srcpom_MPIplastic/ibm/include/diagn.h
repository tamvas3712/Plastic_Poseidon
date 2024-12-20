!--diagnostics
    real ::  tbvariant,tnvariant,tsivariant,weivariant,sizvariant,tbeach,depvariant,denvariant
    real ::  tbvariantBea,tbvariant100,tbvariantBot,tnvariant100,tsivariantBea,tsivariantBot
    real ::  tnvariantBea,tnvariantBot,tsivariant100,fdepvariant,fdepvariantBot,tbvariantBur

    common /blkfish/ tbvariant(nvariantsmax,nspecies)           !total biomass (tn)
    common /blkfish/ tbvariant100(nvariantsmax,nspecies)        !total biomass above 100m depth (tn)
    common /blkfish/ tbvariantBea(nvariantsmax,nspecies)        !total biomass on beach (tn)
    common /blkfish/ tbvariantBur(nvariantsmax,nspecies)        !total biomass buried on beach (tn)
    common /blkfish/ tbvariantBot(nvariantsmax,nspecies)        !total biomass on sea floor (tn)

    common /blkfish/ tnvariant100(nvariantsmax,nspecies)        !mean conc. above 100m depth (#/Km2)
    common /blkfish/ tnvariant(nvariantsmax,nspecies)           !mean conc. above water column (#/Km2)
    common /blkfish/ tnvariantBea(nvariantsmax,nspecies)        !mean conc. on beach (#/m)
    common /blkfish/ tnvariantBot(nvariantsmax,nspecies)        !mean conc  on sea floor (#/Km2)

    common /blkfish/ tsivariant(nvariantsmax,nspecies)          !SIs
    common /blkfish/ tsivariant100(nvariantsmax,nspecies)       !SIs above 100m depth
    common /blkfish/ tsivariantBea(nvariantsmax,nspecies)       !SIs at beach
    common /blkfish/ tsivariantBot(nvariantsmax,nspecies)       !SIs at bottom

    common /blkfish/ fdepvariant(nvariantsmax,nspecies)         !mean depth of sea particles
    common /blkfish/ fdepvariantBot(nvariantsmax,nspecies)      !mean depth of sea floor particles

    common /blkfish/ weivariant(nvariantsmax,nspecies)          !mean weight
    common /blkfish/ sizvariant(nvariantsmax,nspecies)           !mean size
    common /blkfish/ tbeach(nvariantsmax,nspecies)          !total beached
    common /blkfish/ depvariant(nvariantsmax,nspecies)          !mean depth
    common /blkfish/ denvariant(nvariantsmax,nspecies)          !mean density


