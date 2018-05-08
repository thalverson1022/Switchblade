CC=mpif90
FILENAME=SwitchBlade_II_MAIN_v5_1
MODS1=SwitchBlade_II_MODS_v5_1
OUTPUT=switchblade
TAG=f90
LIBS = -Wl,-rpath, ${TACC_MKL_LIB} -L${TACC_SCALAPACK_LIB} -lscalapack ${TACC_SCALAPACK_LIB}/blacs_MPI-LINUX-0.a ${TACC_SCALAPACK_LIB}/blacsF77init_MPI-LINUX-0.a ${TACC_SCALAPACK_LIB}/blacs_MPI-LINUX-0.a -L${TACC_MKL_LIB} -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread 
 

all:execute

execute:$(FILENAME).$(TAG) Makefile
	$(CC) $(MODS1).$(TAG) $(FILENAME).$(TAG) -o $(OUTPUT).exe $(LIBS)
clean:
	rm -f $(OUTPUT).exe; 
	rm -f *.mod;

