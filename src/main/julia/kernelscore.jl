using HPAT

@acc hpat function kernelscore_test(file_name)
    df = DataSource(DataTable{:id=Int64,:x=Float64,:y=Float64}, HDF5, file_name)
    points = [-1.0, 2.0, 5.0]
    N = size(points,1)
    b = 0.5
    exps::Float64 = 0.0 
    @par exps(+) for i in 1:length(df[:x])
        d = -(df[:x][i]-points).^2./(2*b^2)
        m = minimum(d)
        exps += m-log(b*N)+log(sum(exp(d-m)))
    end 
    return exps
end

println(kernelscore_test(ENV["HOME"]*"/tmp/cumsum.hdf5"))

