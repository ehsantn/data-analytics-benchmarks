using HPAT 
using MPI
using DocOpt
ParallelAccelerator.ParallelIR.PIRHoistParfors(1)
ParallelAccelerator.CGen.setRecreateLoops(true)
ParallelAccelerator.CGen.setRecreateConds(true)
HPAT.DistributedPass.DistOptimize(true)
HPAT.DistributedPass.set_debug_level(3)

@acc hpat function logistic_regression(iterations, file_name)
    points = DataSource(Matrix{Float64},HDF5,"/points", file_name)
    responses = DataSource(Vector{Float64},HDF5,"/responses", file_name)
    D,N = size(points) # number of features,samples
    labels = reshape(responses,1,N)
    w = ones(1,D)-.5
    tic()
    for i in 1:iterations
       w -= ((1.0./(1.0.+exp(-labels.*(w*points))).-1.0).*labels)*points'
    end
    toc()
    w
end

function main()
    doc = """Logistic regression statistical method.

Usage:
  logistic_regression.jl -h | --help
  logistic_regression.jl [--iterations=<iterations>] [--file=<file>]

Options:
  -h --help                  Show this screen.
  --iterations=<iterations>  Specify number of iterations; defaults to 20.
  --file=<file>              Specify input file; defaults to HPAT's default generated data file.

"""
    arguments = docopt(doc)
    iterations = 20
    file_name = HPAT.getDefaultDataPath()*"logistic_regression.hdf5"

    if (arguments["--iterations"] != nothing)
        iterations = parse(Int, arguments["--iterations"])
    end

    if (arguments["--file"] != nothing)
        file_name = arguments["--file"]
    end 

    W = logistic_regression(iterations, file_name)

    if MPI.Comm_rank(MPI.COMM_WORLD)==0 println("result = ", W) end
end

main()

