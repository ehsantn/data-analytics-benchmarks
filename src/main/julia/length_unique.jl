using HPAT

@acc hpat function unique_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Int64}, HDF5, file_name)
    df1 = aggregate(df, :id, :sx=length(unique(:x)))
    return sum(df1[:sx])
end

unique_test(ENV["HOME"]*"/tmp/unique.hdf5")

