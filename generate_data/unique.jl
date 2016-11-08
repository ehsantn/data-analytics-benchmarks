using HDF5

function generate_files(path,N)
    hdf5_file = path*"unique.hdf5"
    csv_file = path*"unique.csv"
    if isfile(hdf5_file)
        rm(hdf5_file)
    end
    if isfile(csv_file)
        rm(csv_file)
    end
    # uniform ids in a range for load balance among nodes
    id = rand(1:256, N)
    x = rand(1:256*1000,N)
    h5write(hdf5_file,"/id", id)
    h5write(hdf5_file,"/x", x)
    writecsv(csv_file, zip(id,x))
end

N = 2*10^6
path = ENV["HOME"]*"/tmp/"
generate_files(path, N)
