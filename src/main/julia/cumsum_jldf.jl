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
    tic()
    b = cumsum(df[:x])
    t = toq()
    println("DataFrames.jl filter took: ", t," ", sum(b))
end

main()
