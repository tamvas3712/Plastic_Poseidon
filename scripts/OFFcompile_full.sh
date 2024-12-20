set -x

mkdir ../bin/
rm ../bin/pom.exe
rm ../bin/pomibm.exe
rm ../bin/prepare_netcdfPOM*
rm ../bin/write_netcdfPOM*
rm ../bin/pre_meteo_posII

cd ../srcpom_MPIeco/ecology
rm libecology.a
make clean -f Makefile_med
make -f Makefile_med
cd ..
make clean -f Makefile_eco
make -f Makefile_eco
make clean -f Makefile_eco

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
#compile the full prepN* first, then the JustPlastics
cd ../OLDprepNetcdfplastic/
./compilePlastic
cd ../prepNetcdfplastic/
./compilePlastic
./ifc_ibm
cd ..
echo "done compiling"

cp srcpom_MPIeco/pom.nml bin/
cp srcpom_MPIplastic/pomibm.nml bin/
ln -s srcpom_MPIplastic/ibm/include/ibm.dat data/

echo "done copying settings files and linking"

chmod +775 *  #to be able to work in a cross-user accessible build, comment out if otherwise

cd scripts/

