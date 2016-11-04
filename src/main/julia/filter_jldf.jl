using DataFrames

function main()
    table_path = ARGS[1]
    df = readtable(table_path)
    tic()
    df[df[1].<100,:]
    t = toq()
    println("DataFrames.jl filter took: ", t)
end

main()
