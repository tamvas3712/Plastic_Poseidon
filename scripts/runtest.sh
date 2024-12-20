yr1=`cut -c5-6 dir1`
mo1=`cut -c3-4 dir1`
day1=`cut -c1-2 dir1`

#./../bin/prepare_netcdfPOMwave>../data/wave.out  << EOF
# $day1 $mo1 $yr1
#EOF

./../bin/pre_meteo_posII>../data/out_meteo << EOF
 $day1 $mo1 $yr1
EOF

#./../bin/prepare_netcdfPOMatmos>../data/outmeteo << EOF
# $day1 $mo1 $yr1
#EOF

#time ./../bin/prepare_netcdfPOMibmplasticaddBeachSize>../data/test.out << EOF
# $day1 $mo1
#EOF

#./../bin/prepare_netcdfPOMibmplasticaddWWT>../data/addWWTnew.out << EOF
# $day1 $mo1
#EOF

#./../bin/prepare_netcdfPOMibmplasticaddSou>../data/add.out << EOF
# $day1 $mo1 $yr1
#EOF


#./../bin/write_netcdfPOMibmplastic>../data/plastic.out << EOF
# $day1 $mo1 $yr1
#EOF

#./../bin/prepare_netcdfPOMwave>../data/wave.out  << EOF
# $day1 $mo1 $yr1
#EOF

