#include "par_write.h"

void par_write_wrapper(MPI_Comm comm, MPI_Info info, int mode, std::filesystem::path&& path, size_t dsetsize){
    int rank;
    MPI_Comm_rank(comm, &rank);
    switch(mode){
        case wrt_mode::WRT_HDF5:
            if(rank == 0) std::cout << "Writing using HDF5 APIs" << std::endl;
            par_write_hdf5(comm, info, std::move(path), dsetsize);
            break;
        case wrt_mode::WRT_MPIO:
            if(rank == 0) std::cout << "Writing using MPIIO APIs" << std::endl;
            par_write_mpio(comm, info, std::move(path), dsetsize);
            break;
        case wrt_mode::WRT_POSIX:
            if(rank == 0) std::cout << "Writing using POSIX APIs" << std::endl;
            par_write_posix(comm, info, std::move(path), dsetsize);
            break;
        default:
            if(rank == 0) std::cerr << "Unknown write mode" << std::endl;
            break;
    }
}

void par_write_hdf5(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize){
    using namespace std::chrono;
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
    if (rank == 0)
        std::cout << "Dataset size: " << double(dsetsize * sizeof(int)) / (1024.0 * 1024.0) << "MB" << std::endl;
    size_t start[1]  = {rank * (dsetsize / nprocs)},
           count[1]  = {dsetsize / nprocs},
           stride[1] = {1};
    int* data = new int[*count];
    for (size_t i = 0; i < count[0]; i++)
        data[i] = rank;
    auto h5filename = path / "parwrite.h5";

    ////////////// begin //////////////
    auto begin = system_clock::now();
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    herr_t err;
    H5Pset_fapl_mpio(plist_id, comm, info);

    hid_t fid = H5Fcreate(h5filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    err = H5Pclose(plist_id);
    hid_t sid1 = H5Screate_simple(1, &dsetsize, NULL);
    hid_t dataset = H5Dcreate2(fid, "VAR", H5T_NATIVE_INT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t file_dataspace = H5Dget_space(dataset);
    err = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
    hid_t mem_dataspace = H5Screate_simple(1, count, NULL);

    plist_id = H5Pcreate (H5P_DATASET_XFER); 
    H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_COLLECTIVE);
    err = H5Dwrite(dataset, H5T_NATIVE_INT, mem_dataspace, file_dataspace, plist_id, data);

    H5Sclose(file_dataspace);
    H5Sclose(mem_dataspace);
    H5Pclose(plist_id);
    // err = H5Fflush(dataset, H5F_SCOPE_LOCAL);
    MPI_Barrier(comm);
    /////////////// end ///////////////

    delete[] data;
    auto end = system_clock::now();
    if (rank == 0)
        std::cout << "Parallel HDF5 Write Time cost: " << (double)duration_cast<microseconds>(end - begin).count() / (double)1000000.0 << "s" << std::endl;
}

void par_write_mpio(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize){
    using namespace std::chrono;
    int rank, nprocs;
    auto mpi_filename = path / "mpio.bin";
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
    if (rank == 0)
        std::cout << "Dataset size: " << double(dsetsize * sizeof(int)) / (1024.0 * 1024.0) << "MB" << std::endl;
    if (std::filesystem::exists(mpi_filename))
        std::filesystem::remove(mpi_filename);
    size_t count = dsetsize / nprocs;
    int* data = new int[count];
    for (size_t i = 0; i < count; i++)
        data[i] = rank;
    
    ////////////// begin //////////////
    auto begin = system_clock::now();
    MPI_File fh;
    MPI_Status status;
    int err = MPI_File_open(comm, mpi_filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &fh);
    MPI_File_write_at_all(fh, rank*count*sizeof(int), data, count, MPI_INT, &status);
    MPI_File_close(&fh);
    MPI_Barrier(comm);
    auto end = system_clock::now();
    /////////////// end ///////////////

    delete[] data;
    if (rank == 0)
        std::cout << "MPI-IO Write Time cost: " << (double)duration_cast<microseconds>(end - begin).count() / (double)1000000.0 << "s" << std::endl;
}

void par_write_posix(MPI_Comm comm, MPI_Info info, std::filesystem::path&& path, size_t dsetsize){
    using namespace std::chrono;
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
    if (rank == 0)
        std::cout << "Dataset size: " << double(dsetsize * sizeof(int)) / (1024.0 * 1024.0) << "MB" << std::endl;
    auto begin = system_clock::now();
    

    MPI_Barrier(comm);
    auto end = system_clock::now();
    if (rank == 0)
        std::cout << "POSIX-IO Write Time cost: " << (double)duration_cast<microseconds>(end - begin).count() / (double)1000000.0 << "s" << std::endl;
}