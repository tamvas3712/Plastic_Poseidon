cd ..
rm -r data/Hind*
rm -r data/ECO_RES/*
rm data/log
rm data/logibm

#rm out/*
cp data/ECO_RES_original/ibmplastic010110.nc in/pom.ibmplastic.nc
cp data/ECO_RES_original/restart010110.nc in/restart.nc
cp data/FIELDS.INPUT bin/
cd scripts/
echo "done"

