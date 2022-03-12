#include "posix_empty_grp.h"

static void recursive_mkdir(std::filesystem::path& path, unsigned depth, unsigned nbranches);

void make_posix_empty_grp(std::filesystem::path&& path, unsigned depth, unsigned nbranches){
    using namespace std::chrono;
    auto root_grp = path / "rootgrp/";
    if (depth >= 100){
        std::cerr << "Too deep: " << depth << std::endl;
        return;
    }

    auto total_grps = ((int)nbranches * (1 - std::pow(nbranches, depth))) / (1 - (int)nbranches);
    std::cout << "Total: " << total_grps << " groups" << std::endl;

    if (std::filesystem::exists(root_grp)){
        std::cerr << "Root directory exists, removing..." << std::endl;
        std::filesystem::remove_all(root_grp);
    }
    auto begin = system_clock::now();
    mkdir(root_grp.c_str(), 0755);
    recursive_mkdir(root_grp, depth, nbranches);
    auto end = system_clock::now();
    std::cout << "POSIX Time cost: " << (double)duration_cast<microseconds>(end-begin).count() / (double)1000000.0 << "s" << std::endl;
}

void make_posix_empty_grp_par(std::filesystem::path&& path, unsigned depth, unsigned nbranches){
using namespace std::chrono;
    auto root_grp = path / "rootgrp/";
    int rank, nprocs;
    if (depth >= 100){
        if (rank == 0)
            std::cerr << "Too deep: " << depth << std::endl;
        return;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert(nprocs == nbranches);

    if (rank == 0){
        auto total_grps = ((int)nbranches * (1 - std::pow(nbranches, depth))) / (1 - (int)nbranches);
        std::cout << "Total: " << total_grps << " groups" << std::endl;        
        if (std::filesystem::exists(root_grp)){
            std::cerr << "Root directory exists, removing..." << std::endl;
            std::filesystem::remove_all(root_grp);
        }
    }

    auto begin = system_clock::now();
    if (rank == 0)
        mkdir(root_grp.c_str(), 0755);
    MPI_Barrier(MPI_COMM_WORLD);
    auto proc_path = root_grp / std::to_string(rank);
    mkdir(proc_path.c_str(), 0755);
    recursive_mkdir(proc_path, depth - 1, nbranches);
    MPI_Barrier(MPI_COMM_WORLD);
    auto end = system_clock::now();
    if (rank == 0)
        std::cout << "POSIX Time cost: " << (double)duration_cast<microseconds>(end-begin).count() / (double)1000000.0 << "s" << std::endl;
}

static void recursive_mkdir(std::filesystem::path& path, unsigned depth, unsigned nbranches){
    if (depth == 0)
        return;
    for (auto i = 0; i < nbranches; i++){
        auto subdir = path / std::to_string(i);
        mkdir(subdir.c_str(), 0755);
        recursive_mkdir(subdir, depth - 1, nbranches);
    }
}