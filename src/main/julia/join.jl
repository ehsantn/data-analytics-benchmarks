using HPAT

@acc hpat function cumsum_test(file_name, file_name2)
    df1 = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    df2 = DataSource(DataTable{:id=Int64,:x1=Float64,:y1=Float64}, HDF5, file_name2)
    df3 = join(df1, df2, :id==:id, :id)
    return sum(df3[:x])
end

cumsum_test(ENV["HOME"]*"/tmp/udf.hdf5", ENV["HOME"]*"/tmp/udf2.hdf5")

