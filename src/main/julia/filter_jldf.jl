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
    df[df[:id].<100,:]
    t = toq()
    println("DataFrames.jl filter took: ", t)
end

main()
