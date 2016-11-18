using DataFrames
using HDF5

function main()
    table_path = ARGS[1]
    id = h5read(table_path,"id")
    x = h5read(table_path,"x")
    y = h5read(table_path,"y")

    df = DataFrame()
    df[:id] = id
    df[:x] = x
    df[:y] = y
    A = Array(Float64,size(x,1))
    tic()
    for i in 2:size(x,1)-1
        A[i] = (df[:x][i-1]+df[:x][i]+df[:x][i+1])/3.
    end
    t = toq()
    println("DataFrames.jl SMA took: ", t, " val", sum(A))
end

main()
