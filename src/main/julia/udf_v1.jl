using HPAT

@acc hpat function udf_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    df2 = aggregate(df, :id, :sx = sum(2*:x), :sy = sum(2*:y))
    return sum(df2[:sx])
end

udf_test(ENV["HOME"]*"/tmp/udf.hdf5")

