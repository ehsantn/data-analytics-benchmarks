using HDF5
using DocOpt
using HPAT

@acc hpat function generate_file(hdf5_file,N)
    D = 10
    A = rand(Float64,D,N)
    Y = rand(Float64,N)
    DataSink(A, HDF5,"/points", hdf5_file)
    DataSink(Y, HDF5,"/responses", hdf5_file)
end


function main()
    doc = """generate_logistic_regression.jl

generate random input for logistic regression example.

Usage:
  generate_logistic_regression.jl -h | --help
  generate_logistic_regression.jl [--instances=<instances>] [--path=<path>]

Options:
  -h --help                  Show this screen.
  --instances=<instances>    Specify number of instances; defaults to 2000000.
  --path=<path>              Specify output path for generated files; defaults to HPAT's default data path.
"""

    arguments = docopt(doc)

    if (arguments["--instances"] != nothing)
        instances = parse(Int, arguments["--instances"])
    else
        instances = 2*10^3
    end

    if (arguments["--path"] != nothing)
        path = arguments["--path"]
    else
        path = HPAT.getDefaultDataPath()
    end
    hdf5_file = path*"logistic_regression.hdf5"
    if isfile(hdf5_file)
        rm(hdf5_file)
    end

    generate_file(hdf5_file, instances)
end

main()
