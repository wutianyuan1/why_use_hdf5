#ifndef PAR_WRITE_H
#define PAR_WRITE_H
#include <iostream>
#include <chrono>
#include <filesystem>
#include <hdf5.h>
#include <mpi.h>

enum wrt_mode{ WRT_HDF5=0, WRT_MPIO, WRT_POSIX };

void par_write_wrapper(MPI_Comm comm, MPI_Info info, int mode, std::filesystem::path&& path, size_t dsetsize);

void par_write_hdf5(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize);

void par_write_mpio(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize);

void par_write_posix(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize);

#endif