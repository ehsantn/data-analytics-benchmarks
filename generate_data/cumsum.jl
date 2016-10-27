using HDF5

function generate_files(path,N)
    hdf5_file = path*"cumsum.hdf5"
    csv_file = path*"cumsum.csv"
    if isfile(hdf5_file)
        rm(hdf5_file)
    end
    if isfile(csv_file)
        rm(csv_file)
    end
    # uniform ids in a range for load balance among nodes
    id = collect(1:N)
    x = rand(N)
    y = rand(N)
    h5write(hdf5_file,"/id", id)
    h5write(hdf5_file,"/x", x)
    h5write(hdf5_file,"/y", y)
    writecsv(csv_file, zip(id,x,y))
end

#N = 256*10^6
N = 256
path = ENV["HOME"]*"/tmp/"
generate_files(path, N)
