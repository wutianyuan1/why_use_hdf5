#ifndef WRT_INDEP_VARS_H
#define WRT_INDEP_VARS_H
#include <iostream>
#include <chrono>
#include <filesystem>
#include <hdf5.h>
#include <mpi.h>

enum wrt2_mode{ INDEP_HDF5=0, INDEP_MPIO, INDEP_POSIX };

void indep_write_wrapper(MPI_Comm comm, MPI_Info info, int mode, std::filesystem::path&& path, size_t dsetsize);

void indep_write_hdf5(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize);

void indep_write_mpio(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize);

void indep_write_posix(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize);

#endif //