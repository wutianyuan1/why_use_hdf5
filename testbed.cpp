#include <iostream>
#include <filesystem>
#include <mpi.h>
#include "posix_empty_grp.h"
#include "hdf5_empty_grp.h"
#include "par_write.h"
#include "wrt_indep_vars.h"

void test_empty_grp(int argc, char** argv){
    if (argc <= 3){
        std::cerr << "Usage: ./test <PATH> <DEPTH> <NBRANCHES>" << std::endl;
        exit(1);
    }
    make_hdf5_empty_grp(std::filesystem::path(argv[1]), (unsigned)std::atoi(argv[2]), (unsigned)std::atoi(argv[3]));
    make_posix_empty_grp(std::filesystem::path(argv[1]), (unsigned)std::atoi(argv[2]), (unsigned)std::atoi(argv[3]));
}

void test_empty_grp_par(int argc, char** argv){
    if (argc <= 3){
        std::cerr << "Usage: ./test <PATH> <DEPTH> <NBRANCHES>" << std::endl;
        exit(1);
    }
    make_hdf5_empty_grp_par(std::filesystem::path(argv[1]), (unsigned)std::atoi(argv[2]), (unsigned)std::atoi(argv[3]));
    make_posix_empty_grp_par(std::filesystem::path(argv[1]), (unsigned)std::atoi(argv[2]), (unsigned)std::atoi(argv[3]));
}

void test_par_write(int argc, char** argv){
    if (argc <= 2){
        std::cerr << "Parwrite Usage: ./test <OP> <PATH> <DATASIZE:MB>" << std::endl;
        exit(1);
    }
    double nelems = 1024.0 * 1024.0 * std::atof(argv[3]) / sizeof(int);
    par_write_wrapper(MPI_COMM_WORLD, MPI_INFO_NULL, std::atoi(argv[1]), std::filesystem::path(argv[2]), (size_t)nelems);
}

void test_indep_write(int argc, char** argv){
    if (argc <= 2){
        std::cerr << "Indepwrite Usage: ./test <OP> <PATH> <DATASIZE:MB>" << std::endl;
        exit(1);
    }
    double nelems = 1024.0 * 1024.0 * std::atof(argv[3]) / sizeof(int);
    indep_write_wrapper(MPI_COMM_WORLD, MPI_INFO_NULL, std::atoi(argv[1]), std::filesystem::path(argv[2]), (size_t)nelems);
}

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    test_empty_grp_par(argc, argv);
    MPI_Finalize();
    return 0;
}
