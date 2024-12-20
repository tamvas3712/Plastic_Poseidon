#!/bin/bash
rundir=$1
restdir=01

if test `cat stop_run` -eq 1;then
  mkdir $rundir`cat dir1`

  yr1=`cut -c5-6 dir1`
  mo1=`cut -c3-4 dir1`
  day1=`cut -c1-2 dir1`

  #read atm. fields
./../bin/pre_meteo_posII>../data/out_meteo << EOF
 $day1 $mo1 $yr1
EOF

  #backup restart files
  if test `cut -c1-2 dir1` -eq $restdir ;then
    cp ../in/restart.nc ../data/ECO_RES/restart${day1}${mo1}${yr1}.nc
  fi

  #prepare atmos netcdf
  ./../bin/prepare_netcdfPOMatmos>../data/outmeteo 

  cp ../data/pom.atmos.nc ../in/

  #run hydrodynamic model
  time mpirun -np 1 ./../bin/pom.exe > ../data/log

  #write hydro output and BC's
./../bin/write_netcdfPOMhydro>../data/hydro.out  << EOF
 $day1 $mo1 $yr1
EOF

  mv ../bin/OUTP* $rundir/`cat dir1`
  mv ../bin/OUTB* $rundir/`cat dir1`

  #update restart
  cp ../out/restart.nc ../in/
  cp ../out/pom.hydroibm.nc ../in/


  cp log $rundir/`cat dir1`

  ./makeavscript.sh $rundir

fi


