using HPAT

@acc hpat function movingavg_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    avg = stencil(x->(x[-1]+2*x[0]+x[1])/4.0, df[:x])
    return sum(avg)
end

movingavg_test(ENV["HOME"]*"/tmp/cumsum.hdf5")

