#
# module load nvhpc
#
#CUDA_HOME=~/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/cuda
HDF5_HOME=/opt/hdf5/HDF5-1.12.0-Linux/HDF_Group/HDF5/1.12.0
#
FC = mpif90

#FLAGS = -O4 -ta:tesla,nofma,maxregcount:112 -Minfo=accel -DUSE_CATALYST -DUSE_CUDA -DUSE_NVTX -Mcuda=cuda11.0
FLAGS = -O4
#FLAGS += -DSCALAR  # compute scalar field  

#CXX=mpic++
#CXXFLAGS = -std=c++11

#======================================================================
# Library
#======================================================================
#LINKS += -lstdc++
#INCS=-I$(CUDA_HOME)/include

#LINKS = -lfftw3 -llapack #-lblas #-lmkl_core -lmkl_intel_lp64 -mkl=parallel
#LINKS += -L$(CUDA_HOME)/lib64 -lnvToolsExt -Mcudalib=cufft,nccl,cutensor

LINKS += -L${HDF5_HOME}/lib -lhdf5 -lhdf5_fortran

INCS+=-I${HDF5_HOME}/include
INCS+=-I${HDF5_HOME}/include/static

#====================================================================== 
# Make program
#======================================================================

PROGRAM = field3d

OBJECTS += field3d.o 

#OBJECTS += densmc.o invtrro.o solroi.o solroj.o # compute scalar field
           
#MODULES = param.o decomp_2d.o utils.o nvtx.o

#=======================================================================
# Linking    
#=======================================================================

#$(PROGRAM) : $(OBJECTS) $(MODULES)
#	$(FC) $(FLAGS) $(OBJECTS) $(MODULES) $(LINKS) -o $@ 
	
$(PROGRAM) : $(OBJECTS)
	$(FC) $(FLAGS) $(OBJECTS) $(LINKS) -o $@

##=======================================================================
## Dependencies
##=======================================================================
#
#param.o: param.F90
#	$(FC) $(FLAGS) -c param.F90
#
#decomp_2d.o: decomp_2d.F90
#	$(FC) $(FLAGS) -c decomp_2d.F90
#
#utils.o: utils.F90
#	$(FC) $(FLAGS) -c utils.F90
#
#nvtx.o: nvtx.F90
#	$(FC) $(FLAGS) -c nvtx.F90
#
%.o:   %.f90 $(MODULES)
	$(FC) $(FLAGS) $(INCS) -c $<

#%.o: %.cxx $(MODULES)
#	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
#
##=======================================================================
## Clean up
##=======================================================================
#
#clean : 
#	rm *.o *.mod 
