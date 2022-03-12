#include "hdf5_empty_grp.h"

static void recursive_crt_grp(hid_t file, unsigned depth, unsigned nbranches);

void make_hdf5_empty_grp(std::filesystem::path&& path, unsigned depth, unsigned nbranches){
    using namespace std::chrono;
    auto h5filename = path / "test.h5";
    if (depth >= 100){
        std::cerr << "Too deep: " << depth << std::endl;
        return;
    }

    auto total_grps = ((int)nbranches * (1 - std::pow(nbranches, depth))) / (1 - (int)nbranches);
    std::cout << "Total: " << total_grps << " groups" << std::endl;

    auto begin = system_clock::now();
    hid_t h5file = H5Fcreate(h5filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    recursive_crt_grp(h5file, depth, nbranches);
    H5Fclose(h5file);
    auto end = system_clock::now();
    std::cout << "HDF5 Time cost: " << (double)duration_cast<microseconds>(end-begin).count() / (double)1000000.0 << "s" << std::endl;
}

void make_hdf5_empty_grp_par(std::filesystem::path&& path, unsigned depth, unsigned nbranches){
    using namespace std::chrono;
    auto h5filename = path / "test.h5";
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (depth >= 100){
        if (rank == 0)
            std::cerr << "Too deep: " << depth << std::endl;
        return;
    }
    
    auto total_grps = ((int)nbranches * (1 - std::pow(nbranches, depth))) / (1 - (int)nbranches);
    if (rank == 0)
        std::cout << "Total: " << total_grps << " groups" << std::endl;

    auto begin = system_clock::now();
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t h5file = H5Fcreate(h5filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    recursive_crt_grp(h5file, depth, nbranches);
    H5Fclose(h5file);
    MPI_Barrier(MPI_COMM_WORLD);
    auto end = system_clock::now();
    if (rank == 0)
        std::cout << "HDF5 Time cost: " << (double)duration_cast<microseconds>(end-begin).count() / (double)1000000.0 << "s" << std::endl;
}

static void recursive_crt_grp(hid_t file, unsigned depth, unsigned nbranches){
    if (depth == 0)
        return;
    for (auto i = 0; i < nbranches; i++){
        auto grpname = std::to_string(i);
        auto gid = H5Gcreate(file, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        recursive_crt_grp(gid, depth - 1, nbranches);
        H5Gclose(gid);
    }
}