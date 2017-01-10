using HPAT 
using MPI
using DocOpt
ParallelAccelerator.ParallelIR.PIRHoistParfors(1)
ParallelAccelerator.CGen.setRecreateLoops(true)
ParallelAccelerator.CGen.setRecreateConds(true)
HPAT.DistributedPass.DistOptimize(true)
#HPAT.DistributedPass.set_debug_level(3)

@acc hpat function logistic_regression(iterations, N)
    D = 10  # Number of features
    labels = rand(1,N)
    points = rand(D,N)
    w = reshape(2.0.*rand(D)-1.0,1,D)
    tic()
    for i in 1:iterations
       w -= ((1.0./(1.0.+exp(-labels.*(w*points))).-1.0).*labels)*points'
    end
    toc()
    w
end

function main()
    doc = """logistic_regression.jl

Logistic regression statistical method.

Usage:
  logistic_regression.jl -h | --help
  logistic_regression.jl [--iterations=<iterations>] [--instances=<instances>]

Options:
  -h --help                  Show this screen.
  --iterations=<iterations>  Specify number of iterations; defaults to 20.
  --instances=<instances>    Specify number of instances; defaults to 10^7.
"""
    arguments = docopt(doc)
    iterations = 20
    instances = 10^7

    if (arguments["--iterations"] != nothing)
        iterations = parse(Int, arguments["--iterations"])
    end

    if (arguments["--instances"] != nothing)
        instances = parse(Int, arguments["--instances"])
    end

    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    W = logistic_regression(iterations, instances)

    if MPI.Comm_rank(MPI.COMM_WORLD)==0 println("result = ", W) end
end

main()

