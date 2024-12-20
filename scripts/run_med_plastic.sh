#!/bin/bash
set -x #comment if you dont want debugging mode where each line of this script is printed during runtime
rundir=$1
restdir=01

if test `cat stop_run` -eq 1;then
  mkdir $rundir`cat dir1`

  yr1=`cut -c5-6 dir1`
  mo1=`cut -c3-4 dir1`
  day1=`cut -c1-2 dir1`



  cd ../bin


  #read atm. fields PLACEHOLDER
#./pre_meteo_posII>../data/out_meteo << EOF
# $day1 $mo1 $yr1
#EOF


  #backup restart files
#  if test `cut -c1-2 dir1` -eq $restdir ;then
#    cp ../in/restart.nc ../data/ECO_RES/restart${day1}${mo1}${yr1}.nc
#    cp ../in/pom.ibmplastic.nc ../data/ECO_RES/ibmplastic${day1}${mo1}${yr1}.nc
#  fi


  #prepare atmos netcdf PLACEHOLDER
#  ./prepare_netcdfPOMatmos>../data/outmeteo 
#  cp pom.atmos.nc ../in/


  #run hydrodynamic model PLACEHOLDER
#  time mpirun -np 3 ./pom.exe > ../data/log


  #write hydro output and BC's PLACEHOLDER
#./write_netcdfPOMhydro>../data/hydro.out  << EOF
# $day1 $mo1 $yr1
#EOF
#  mv OUTP* $rundir/`cat dir1`
#  mv OUTB* $rundir/`cat dir1`



  #update restart  PLACEHOLDER
#  cp ../out/restart.nc ../in/    #normally both are a mv command and not cp
#  cp ../out/pom.hydroibm.nc ../in/



  #prepare wave netcdf PLACEHOLDER
#./prepare_netcdfPOMwave>../data/wave.out  << EOF
  # $day1 $mo1 $yr1
#EOF
#  mv pom.wave.nc ../in/.


  #run plastics model
  time mpirun -np 5 ./pomibm.exe > ../data/logibm


  #write plastics output
./write_netcdfPOMibmplastic>../data/plastic.out << EOF
 $day1 $mo1 $yr1
EOF
  #cp boundarySI* $rundir/`cat dir1`
mv FISH*UTC ${rundir}/FISHOUT/.


  #add plastics sources
./prepare_netcdfPOMibmplasticaddSources>../data/add.out << EOF
   $day1 $mo1
EOF


  #update plastics restart file
  cp ../data/pom.ibmplasticNew.nc ../in/pom.ibmplastic.nc

#cp ../out/pom.ibmplastic.nc ../in/pom.ibmplastic.nc

  cp ../data/log $rundir/`cat dir1`
  cp ../data/logibm $rundir/`cat dir1`
  cp ../data/add*.out $rundir/`cat dir1`
  cp ../data/plastic.out $rundir/`cat dir1`


   # read ecology fields?? PLACEHOLDER
#  cd ../scripts
#  ./makeavscript.sh $rundir

fi

echo "finished with run_med_plastic.sh"
