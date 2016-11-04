using HDF5

function generate_files(path,N)
    hdf5_file = path*"filter_2b.hdf5"
    csv_file = path*"filter_2b.csv"
    if isfile(hdf5_file)
        rm(hdf5_file)
    end
    if isfile(csv_file)
        rm(csv_file)
    end
    # uniform ids in a range for load balance among nodes
    id = rand(1:256, N)
    x = rand(N)
    y = rand(N)
    h5write(hdf5_file,"/id", id)
    h5write(hdf5_file,"/x", x)
    h5write(hdf5_file,"/y", y)
    writecsv(csv_file, zip(id,x,y))
end

N = 2*4*256*10^6
path = ENV["HOME"]*"/tmp/"
generate_files(path, N)
