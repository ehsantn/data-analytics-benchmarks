using HPAT

@acc hpat function cumsum_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    df1 = aggregate(df, :id, :sx=sum(:x), :sy=sum(:y))
    return sum(df1[:sx])
end

cumsum_test(ENV["HOME"]*"/tmp/udf.hdf5")

