using HDF5

function generate_files(path,N)
    hdf5_file = path*"join_500k.hdf5"
    csv_file = path*"join_500k.csv"
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

N = 5*10^5
path = ENV["HOME"]*"/tmp/"
generate_files(path, N)
