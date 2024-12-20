set -x

mkdir ../bin/
rm ../bin/pomibm.exe
rm ../bin/prepare_netcdfPOM*
rm ../bin/write_netcdfPOM*

cd ../srcpom_MPIplastic/ibm/
rm libibm.a
make clean -f Makefile_medplastic
make -f Makefile_medplastic

cd ..
make clean -f Makefile_plastic
make -f Makefile_plastic
make clean -f Makefile_plastic

cd ../prepare/
./ifc_compHyd
./ifc_pre_meteo

cd ../prepNetcdfplastic/
./compilePlastic
./ifc_ibm

cd ..
echo "done compiling"

cp srcpom_MPIplastic/pomibm.nml bin/
ln -s srcpom_MPIplastic/ibm/include/ibm.dat data/

echo "done copying settings files and linking"

chmod +775 *  #to be able to work in a cross-user accessible build, comment out if otherwise

cd scripts/
