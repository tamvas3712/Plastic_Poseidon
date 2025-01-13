# Plastic\_Poseidon 

## Description
Simulator of sea plastic polution, based on the Princeton Ocean Model (POM), ERSEM, and a langrangian particles description of micro and macroplastics. The model is is set up for the Mediterranean sea, but can be adopted for other regions. 

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Authors](#authors)
- [License](#license)
- [General](#general)
- [NetCDF](#netcdf)
- [Plastics](#plastics)
- [Run](#run)

## Installation
 The simulator uses netcdf libraries and the ifort compiler, and normally runs on Intel CPUs. You can use the following commands to install
 dependencies and the simulator. 
 The repository is missing the /data folder, where all netcdf and other data files reside. To get the data folder for the moment please contact us.
```bash
#install make:
sudo apt install make

#install netcdf libraries:
sudo apt-get install netcdf-bin
sudo apt-get install libnetcdf-dev libnetcdff-dev
nc-config --all

#install ifort (see Intel ifort homepage, change filenames accordingly):
sudo sh l_dpcpp-cpp-compiler_p_2022.0.2.84_offline.sh
sudo sh l_fortran-compiler_p_2022.0.2.83_offline.sh

sudo sh ./l_HPCKit_p_2022.1.2.117_offline.sh
sudo sh ./l_BaseKit_p_2022.1.2.146_offline.sh

#change oneapi variables to reflect architecture:
source /opt/intel/oneapi/setvars.sh intel64
sudo apt install build-essential

#(optional) install GrADS visualiser:
sudo apt install grads

#download and compile Plastic_Poseidon:
git clone https://github.com/tamvas3712/Plastic_Poseidon
cd Plastic_Poseidon/scripts/
./compile_everything.sh
```

## Usage
See also [General](#general)
```bash
cd Plastic_Poseidon/scripts/
./reset_simulation.sh
./runscript
```


## Documentation
We use FORD for documentation. Please install [FORD](https://forddocs.readthedocs.io/en/latest/) and generate the documentation using the commands:
```bash
pip install ford
cd Plastic_Poseidon/scripts/
./make_documentation.sh
```

## Authors
George Triantafyllou, Kostas Tsiaras, Annika Polani, Ioannis Tamvakis, others (HCMR).

## License
This project should only be used with the consent of George Triantafyllou and Kostas Tsiaras (HCMR).


## General
Please see the scripts mentioned before to get an idea of how the simulator is compiled and run.

The main script is runscript.sh where we have nested loops to run the simulation for every day of every month of every year. It then calls, in this particular version of the simulator, the run\_med\_plastic.sh (hydrodynamic and ibm model) for every simulated day. Stores results in datefolders in data/.
Specify start/end date and directory to save output (rundir) 

In run\_med\_plastic.sh then, if the file stop\_run has 1 in it, 
  will load atmospheric fields using ./pre\_meteo\_posII, compiled using 
```bash
ifort -o pre_meteo_posII pre_meteo_posII.f90
```
  backup restart files,
  prepare atmos netcdf ./bin/prepare\_netcdfPOMatmos>../data/outmeteo,
  run hydrodynamics
  write hydro out
 
In makeavscript.sh you specify whether you want a 10-day, monthly etc average
  this runs an executable (e.g. make\_hydav) that is created in prepare/ (see prepare/ifc\_comp)

 
To see the netcdf files use the command (example):
```bash
  ncdump -h in/pom.grid.nc 
```

TODO:
All .h files are normally linked from src\*/area/ directories, to be used accordingly for the compilation of codes in other directories, for example prepN\*/. This is destroyed by the copying of the files from the cluster. I need to fix this so that we dont have too many of the same file. Also, all files that are used by "include" statements do not exist in the same folder so FORD documenter does not parse the \*90 files correctly.




## NetCDF
prepare netCDF files to be used by parallel code.

For compilation execute:
```bash
cd prepNetcdfplastic/ 
./compilePlastic
```

a) Pre-processing  
prepare\_netcdfPOMibmplastic.f90 -> prepare initial plastics distribution(uniform), reads data/INITnew.DAT binary (model grid)
Then copy the resulting pom.ibmplastic.nc in directory in/.

b) Runtime  
prepare\_netcdfPOMatmos.f90  -> this reads data/FIELDS.INPUT that is created (each day) by pre\_meteo\_posII.f90 (see scripts/run\_med\_) and creates the netcdf file pom.atmos.nc
prepare\_netcdfPOMwave.f90   -> reads wave output and creates  the netcdf file pom.wave.nc
write\_netcdfPOMhydro.f90 -> reads model output from out/pom.hydroibm.nc and write daily binary files and BCs for nested domains that are stored in rundir (specified in runscript)
write\_netcdfPOMibmplastic.f90 -> reads model output from out/pom.ibmplastic.nc and write daily binary files (FISH\*) that are stored in rundir/FISHOUT/




## Parallel
Parallel domain-decomposition code 

The main idea for domain-decomposition is to use netcdf for all input files(forcing, initial conditions etc)
so that each CPU can read its own domain. During runtime there is exchange of information at the boundaries of
the neiboughring (overlapping) domains 
e.g if domain1 extends from i=1-11 and domain2 from i=10-20 of the (global) model grid, 
then the solution at i=10 is obtained from domain1 and the solution at i=11 from domain2. 
Thus, for some calculated array in the two domains (X1, X2) (both having i=1-11 points on their local grid)
we have X2(1)=X1(10) & X1(11)=X2(2) (see code in pom/parallel\_mpi.f)

To get identical results between the parallel and sequential codes one has to exchange all necessary arrays
at the boundaries of domains. 

Steps to configure/run the model:  
1) Prepare necessary files in netcdf (model grid, ecology initial conditions, hydrodynamic fields etc)  
   see prepNetcdf/README (you also need to create some symbolic links of executables in the "running directory)  

2) Specify the number of domains/cpus in srcpom\_MPI\*/pom.h  
   e.g if The global domain is  
     $  im\_global=872  ,  
     $  jm\_global=317  ,  
     $  kb=25          ,  

 You specify the part of the grid each cpu calculates, i.e. you calculate im\_local,jm\_local  
 To do this you take (im\_global-2,jm\_global-2) and divide by a number (n\_procX, n\_procY) that gives an integer  
 Example:  
     $  im\_local=292    ,  
     $  jm\_local=65    ,  
     $  n\_proc=15       )  

  im\_local=(im\_global-2)/3+2=292  
  im\_local=(im\_global-2)/5+2=65  
  n\_proc=n\_procX\*n\_procY=3\*5=15  


  For a different application of the code (e.g different domain etc) main changes are found in   
  srcpom\_MPI\*/area/med20/\*.h. The rest of the code can be more or less the same.


3)Compile the code (see scripts/compile\_everything.sh)

  a) hydrodynamics (please provide your own, or contact us for this part of the codebase)
```bash
cd srcpom_MPIeco
make -f Makefile_eco
```
  This will create executable bin/pom.exe

  b) ibm (in this case a chunk of the total number of SIs (lagrangian particles) is assigned to the specified #CPUs
```bash
cd srcpom_MPIplastic/ibm/
make -f Makefile_medplastic clean
make -f Makefile_medplastic
```
  (when you change the number of domains you need to recompile in srcpom\_MPIibm\*/ibm/)

```bash
cd ..
make -f Makefile\_eco clean
make -f Makefile_eco
cd ..
```
This will create executable bin/pomibm.exe




## Plastic
program main: srcpom\_MPIplastic/pom/pom.f90 (calls initialize.F90 & advance.F90)

Plastics main code: srcpom\_MPIplastic/ibm/lagr\_move.F90
   Uses routines IBM\*.F90 (convert coordinates, beaching, interpolation etc)

The plastics model runs along-side a hydrodynamic model (srcpom\_MPIhydro) providing daily 
(every 3-hour) fields (velocity, vertical/horizontal diffusion etc). 

 
