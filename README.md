## MAGIC master tutorial

This version of MAGIC has the capability to run with both time-dependent profiles and static profiles. The majority of the code remains the same between the two versions, however, the time dependent profile version utilizes setaux_dynamicbg.f and /mpi/claw3ez_mpi_dynamicbg.f, whereas the static version uses setaux_staticbg.f and /mpi/claw3ez_mpi_staticbg.f. There are two version of the MakeFile that can be called to compile either version (if using the non CMake method to build).

### To compile with MakeFile:
Make xclawmpihdf
For the static version (uses just one profile.data file)

Make -f Makefile.dynamic xclawmpihdf
For the dynamic version (uses a profile9XXX.data file for each time step, i.e. first time step is profile9000.data, then profile9001,data etc).

This version also has both the thunderstorm forcing for precipitation driven files and the analytical plume forcing. To run a precipitation driven case uses the flag forcemth =2 in setprob.data. For an analytical plume use forcemth=3.

### To compile with CMAKE:

Firstly, you need h5fortran installed, this can be done using the following:
-Cd into the h5fortran/scripts directory (you may also need to load cmake first)

cmake -B build -DCMAKE_INSTALL_PREFIX=~/local -Dhdf5_parallel=on

cmake --build build

Once this is installed, go into the MAGIC_master directory and type:

cmake -B build -DHDF5_ROOT=~/local

cmake --build build

This will build the version with the static MAGIC backgrounds. To build the dynamic background version, cp CMakeLists_dynamic.txt to CMakeLists.txt, remove the build folder ( rm -rf build) and follw the steps above.

### Example Cases:

#### 1.	Analytically specified plume

The analytical specified plume example is in the Analytical_Plume_example folder. To run this case, first compile the main code in the MAGIC_master directory using either the makefile or CMAKE methods described above, then copy the resulting xclawmpihdf executable to the Analytical_Plume_ example folder.
Next, run the build_ambient_profile.m script in the Matlab_case_set_up folder which will produce the ambient profile file named profile.data (options for the profile can be changed in the build_ambient_profile.m script). Move this profile.data file to the Analytical_Plume_example directory.

The case options are set in claw3ez.data and setprob.data. Claw3ez.data sets up the domain and time while setprob.data sets up the source. Both of these can be edited to suit your needs, but the default is a 500km x 400km x 300km with 2 x 2 x 1 km resolution. The simulation runs for 6000 seconds with outputs every minute. The source has a gaussian half width of 6km in x and y, and 600 secs in time (centered at x=250km, y=200km, and t=2400secs) and has a peak amplitude of 20 mm/ 10 minutes precipitation rate.

To run the cases on HPC queues, use:
qsub pbs_script_DOD if on Narwhal DoD HPC

or msub pbs_script_VEGA if on VEGA.

Ultimately, the command to run the case is:
mpirun -np 125 ./xclawmpi

#### 2.	Real radar driven case
