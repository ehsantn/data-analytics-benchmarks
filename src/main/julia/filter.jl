using HPAT

@acc hpat function cumsum_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    df1 = df[:id<100]
    return sum(df1[:x])
end

cumsum_test(ENV["HOME"]*"/tmp/udf.hdf5")

