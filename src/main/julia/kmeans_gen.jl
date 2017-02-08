using HPAT 
using MPI
using DocOpt
using ParallelAccelerator
ParallelAccelerator.ParallelIR.PIRHoistParfors(1)
#ParallelAccelerator.ParallelIR.set_debug_level(3)
ParallelAccelerator.CGen.setRecreateLoops(true)
ParallelAccelerator.CGen.setRecreateConds(true)
HPAT.DistributedPass.DistOptimize(true)
#HPAT.DistributedPass.set_debug_level(3)

@acc hpat function kmeans(numCenter, iterNum, N)
    D = 10  # Number of features
    points = rand(D,N) 
    centroids = rand(D, numCenter)

    for l in 1:iterNum
        dist::Vector{Vector{Float64}} = [ Float64[sqrt(sum((points[:,i]-centroids[:,j]).^2)) for j in 1:numCenter] for i in 1:N]
        labels::Vector{Int} = [indmin(dist[i]) for i in 1:N]
        centroids::Matrix{Float64} = [ sum(points[j,labels.==i])/sum(labels.==i) for j in 1:D, i in 1:numCenter]
    end 
    return centroids
end

function main()
    doc = """K-means clustering algorithm.

Usage:
  kmeans.jl -h | --help
  kmeans.jl [--iterations=<iterations>] [--instances=<instances>] [--centers=<centers>]

Options:
  -h --help                  Show this screen.
  --iterations=<iterations>  Specify number of iterations; defaults to 20.
  --instances=<instances>    Specify number of instances; defaults to 10^7.
  --centers=<centers>        Specify number of centers; defaults to 5.
"""
    arguments = docopt(doc)

    iterations = 20
    instances = 10^7
    numCenter = 5
    
    if (arguments["--iterations"] != nothing)
        iterations = parse(Int, arguments["--iterations"])
    end

    if (arguments["--instances"] != nothing)
        instances = parse(Int, arguments["--instances"])
    end

    if (arguments["--centers"] != nothing)
        numCenter = parse(Int, arguments["--centers"])
    end

    centroids_out = kmeans(numCenter, iterations, instances)

    if MPI.Comm_rank(MPI.COMM_WORLD)==0 println("result = ", centroids_out) end
end

main()
