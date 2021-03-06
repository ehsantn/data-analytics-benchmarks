using HDF5
using DocOpt
using HPAT
#HPAT.DistributedPass.set_debug_level(3)
#HPAT.DomainPass.set_debug_level(3)
#HPAT.CaptureAPI.set_debug_level(3)
#HPAT.DomainPass.set_debug_level(3)
#using ParallelAccelerator
#ParallelAccelerator.CGen.set_debug_level(3)
#ParallelAccelerator.DomainIR.set_debug_level(3)

@acc hpat function generate_file(hdf5_file, txt_file, N)
    D = 10
    A = randn(D,N)
    w = rand(1,D)
    print(w)
    S = (w*A .> .0)+0.0
    B = ParallelAccelerator.API.reshape(S, (N,))
    T = rand(N) .> .95
    Y = map((t,b)-> t?1.0-b:b, T, B)
    #Y = ones(N).+10.0
    #@par for i in 1:N
    #    Y[i] = T[i] ? B[i] : 1-B[i]
    #end
    DataSink(A, HDF5,"/points", hdf5_file)
    DataSink(Y, HDF5,"/responses", hdf5_file)
    points = [A; ParallelAccelerator.API.reshape(Y,1,N)]
    DataSink(points, TXT, txt_file)
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
    txt_file = path*"logistic_regression.csv"
    if isfile(txt_file)
        rm(txt_file)
    end

    generate_file(hdf5_file, txt_file, instances)
end

main()
