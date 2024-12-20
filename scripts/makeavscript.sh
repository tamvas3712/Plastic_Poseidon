rundir=$1

monthd1=(31 28 31 30 31 30 31 31 30 31 30 31)
monthd2=(31 29 31 30 31 30 31 31 30 31 30 31)


yea=`cut -c5-6 dir1`

  case $yea in
    80) year=1980;;
    81) year=1981;;
    82) year=1982;;
    83) year=1983;;
    84) year=1984;;
    85) year=1985;;
    86) year=1986;;
    87) year=1987;;
    88) year=1988;;
    89) year=1989;;
    90) year=1990;;
    91) year=1991;;
    92) year=1992;;
    93) year=1993;;
    94) year=1994;;
    95) year=1995;;
    96) year=1996;;
    97) year=1997;;
    98) year=1998;;
    99) year=1999;;
    00) year=2000;;
    01) year=2001;;
    02) year=2002;;
    03) year=2003;;
    04) year=2004;;
    05) year=2005;;
    06) year=2006;;
    07) year=2007;;
    08) year=2008;;
    09) year=2009;;
    10) year=2010;;
    11) year=2011;;
    12) year=2012;;
    13) year=2013;;
    14) year=2014;;
    15) year=2015;;
    16) year=2016;;
    17) year=2017;;
    18) year=2018;;
    19) year=2019;;
    20) year=2020;;
    21) year=2021;;
    22) year=2022;;
  esac

echo $year

if test `cut -c3-3 dir1` -eq 0;then
imo=`cut -c4-4 dir1`
else
imo=`cut -c3-4 dir1`
fi


if test `expr $year % 4` -eq 0; then
mdays=${monthd2[${imo}-1]}
else
mdays=${monthd1[${imo}-1]}
fi

outdays=(10 20 $mdays)


if test `cut -c1-1 dir1` -eq 0;then
iday=`cut -c2-2 dir1`
else
iday=`cut -c1-2 dir1`
fi

flag=0
for i in 1 2 3
do 
if test $iday -eq ${outdays[${i}-1]};then
flag=1
fi
done

flag=0
if test $iday -eq $mdays ;then
flag=1
fi


if test $flag -eq 1;then

if test $imo -lt 10;then
imo=0$imo
fi


echo MAKING_AVERAGE
echo $iday
echo $imo
echo $mdays
echo $rundir
./../bin/make_hydav << EOF
 $iday $imo $year $rundir
EOF

if test $iday -eq $mdays ;then
rm ${rundir}*$imo`cut -c5-10 dir1`/EOUT*
rm ${rundir}*$imo`cut -c5-10 dir1`/ECOR*
rm ${rundir}*$imo`cut -c5-10 dir1`/OUTP*
fi



fi

#if test `cut -c1-2 dir1` -eq $restdir ;then
#grep STARTING_DATE ../../MED/inp.prm >test.prm
#echo 'DIRECTORY='`pwd` >>test.prm 

