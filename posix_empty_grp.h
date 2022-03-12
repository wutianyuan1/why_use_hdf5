#ifndef POSIX_EMPTY_GRP_H
#define POSIX_EMPTY_GRP_H
#include <iostream>
#include <filesystem>
#include <string>
#include <chrono>
#include <cmath>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <assert.h>
#include <mpi.h>

void make_posix_empty_grp(std::filesystem::path&& path, unsigned depth, unsigned nbranches);

void make_posix_empty_grp_par(std::filesystem::path&& path, unsigned depth, unsigned nbranches);

#endif //