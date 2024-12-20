#!/bin/bash
#SPECIFY start and end date

year=2010
imoINIT=1
idayINIT=1

yearEND=2010
imoEND=5
idayEND=31 #31

init=1
endyear=0
endmo=0

#SPECIFY name of Run
rundir=../data/HindPlastic31/
mkdir $rundir
mkdir $rundir/FISHOUT/

months=(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
monthd1=(31 28 31 30 31 30 31 31 30 31 30 31)
monthd2=(31 29 31 30 31 30 31 31 30 31 30 31)

while test $year -le $yearEND;
do

  #rm ../data/sst8day.bin
  #ln -s ../prepNetcdfassim/sst8dayMed10_${year}.bin ../data/sst8day.bin

  yyear=`expr $year - 2000`

  if test $yyear -lt 0;then
    yyear=`expr $year - 1900`
  fi

  if test $yyear -lt 10; then
    iyear=0$yyear
  else
    iyear=$yyear
  fi

  if test $init -eq 1; then
    imo=$imoINIT
  else
    imo=1
  fi
  if test $year -eq $yearEND; then
    endyear=1
  fi

  if test $endyear -eq 1; then
    imoyr=$imoEND
  else
    imoyr=12
  fi

  while test $imo -le $imoyr;
  do
    echo $imo
    if test `expr $year % 4` -eq 0; then
      mdays=${monthd2[${imo}-1]}
    else
      mdays=${monthd1[${imo}-1]}
    fi
    echo $mdays
    if test $imo -eq $imoEND -a $endyear -eq 1; then
      endmo=1
    fi

    if test $init -eq 1; then
      iday=$idayINIT
    else
      iday=1
    fi

    init=0

    if test $endmo -eq 1; then
      idaymo=$idayEND
    else
      idaymo=$mdays
    fi

    while test $iday -le $idaymo;
    do
      if test $imo -lt 10; then
        iimo=0$imo
      else
        iimo=$imo
      fi

      if test $iday -lt 10; then
        iiday=0$iday
      else
        iiday=$iday
      fi

      echo ${iiday}${iimo}${iyear}_OCE>dir1
      echo ${iiday}${iimo}${iyear}_OCE>../bin/dir1
      echo ${iiday}${iimo}${iyear}

      month=${months[${imo}-1]}
   
      # executable script
      echo "calling run_med_plastic"
      ./run_med_plastic.sh $rundir

      iday=`expr $iday + 1`
    done
    imo=`expr $imo + 1`
  done
  year=`expr $year + 1`
done

