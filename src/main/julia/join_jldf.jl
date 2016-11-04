using DataFrames
using HDF5

function main()
    table_path = ARGS[1]
    id = h5read(table_path,"id")
    x = h5read(table_path,"x")
    y = h5read(table_path,"y")

    id1 = h5read(table_path,"id")
    x1 = h5read(table_path,"x")
    y1 = h5read(table_path,"y")

    df1 = DataFrame()
    df1[:id] = id
    df1[:x] = x
    df1[:y] = y

    df2 = DataFrame()
    df2[:id] = id1
    df2[:x1] = x1
    df2[:y1] = y1

    tic()
    join(df1,df2,on=:id)
    t = toq()

    println("DataFrames.jl join took: ", t)
end

main()
