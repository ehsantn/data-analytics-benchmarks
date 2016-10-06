using HPAT

HPAT.DistributedPass.set_debug_level(3)

@acc hpat function udf_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    myudf = x->2*x
    df2 = aggregate(df, :id, :sx = sum(myudf(:x)), :sy = sum(myudf(:y)))
    return sum(df2[:sx])
end

udf_test(ENV["HOME"]*"/tmp/udf.hdf5")

