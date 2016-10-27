using HPAT

@acc hpat function cumsum_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    cs = cumsum(df[:x])
    return sum(cs)
end

cumsum_test(ENV["HOME"]*"/tmp/cumsum.hdf5")

