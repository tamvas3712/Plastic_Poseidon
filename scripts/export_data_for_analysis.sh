cd ..
cd data/HindPlastic31/FISHOUT/
ls FISH*>listdays
cd ../../../prepNetcdfplastic/
./ifc_ibm
./readIBM
echo "done, look for popAll*seas*.bin file in prepNetcdf*/"

