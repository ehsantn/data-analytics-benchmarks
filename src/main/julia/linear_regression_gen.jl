using HPAT 
using MPI
using DocOpt
ParallelAccelerator.ParallelIR.PIRHoistParfors(1)
ParallelAccelerator.CGen.setRecreateLoops(true)
ParallelAccelerator.CGen.setRecreateConds(true)
HPAT.DistributedPass.DistOptimize(true)
#HPAT.DistributedPass.set_debug_level(3)

@acc hpat function linear_regression(iterations, N)
    D = 10  # Number of features
    p = 4
    labels = rand(p,N)
    points = rand(D,N)
    w = zeros(p,D)
    alphaN = 0.01/N

    for i in 1:iterations
       w -= alphaN*((w*points)-labels)*points'
    end
    w
end

function main()
    doc = """linear regression statistical method.

Usage:
  linear_regression.jl -h | --help
  linear_regression.jl [--iterations=<iterations>] [--instances=<instances>]

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

    W = linear_regression(iterations, instances)

    if MPI.Comm_rank(MPI.COMM_WORLD)==0 println("result = ", W) end
end

main()

