CXX=mpicxx

# Your HDF5 installation path
HDF5_DIR=/home/wuty/ncinstall

INCLUDES=-I$(HDF5_DIR)/include
LIBS=-L$(HDF5_DIR)/lib -lhdf5 -lstdc++fs
CXXFLAGS=$(INCLUDES) -O3 -std=c++17

TARGET=test

all: testbed

hdf5_grp: hdf5_empty_grp.cpp
	$(CXX) $(CXXFLAGS) -c hdf5_empty_grp.cpp

posix_grp: posix_empty_grp.cpp
	$(CXX) $(CXXFLAGS) -c posix_empty_grp.cpp

par_write: par_write.cpp
	$(CXX) $(CXXFLAGS) -c par_write.cpp

indep_vars: wrt_indep_vars.cpp
	$(CXX) $(CXXFLAGS) -c wrt_indep_vars.cpp

testbed: hdf5_empty_grp.o posix_empty_grp.o par_write.o wrt_indep_vars.o
	$(CXX) $(CXXFLAGS) -o $(TARGET) hdf5_empty_grp.o posix_empty_grp.o par_write.o wrt_indep_vars.o testbed.cpp $(LIBS)

clean:
	rm -f $(TARGET) *.o
	rm -f *.h5
	rm -f *.bin
