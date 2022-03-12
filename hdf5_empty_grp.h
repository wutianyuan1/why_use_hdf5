#ifndef HDF5_EMPTY_GRP_H
#define HDF5_EMPTY_GRP_H
#include <iostream>
#include <filesystem>
#include <string>
#include <chrono>
#include <cmath>
#include <hdf5.h>
#include <mpi.h>

void make_hdf5_empty_grp(std::filesystem::path&& path, unsigned depth, unsigned nbranches);
void make_hdf5_empty_grp_par(std::filesystem::path&& path, unsigned depth, unsigned nbranches);

#endif //