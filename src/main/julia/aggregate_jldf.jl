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
    df2 = aggregate(df, :id, sum)
    t = toq()
    println("DataFrames.jl aggregate took: ", t)
end

main()
