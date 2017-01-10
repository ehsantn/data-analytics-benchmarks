using HPAT
using ParallelAccelerator
ParallelAccelerator.ParallelIR.PIRHoistParfors(1)
ParallelAccelerator.CGen.setRecreateLoops(true)
ParallelAccelerator.CGen.setRecreateConds(true)

@acc hpat function kernelscore_test(n)
    X = rand(n)
    points = [-1.0, 2.0, 5.0]
    N = size(points,1)
    b = 0.5
    exps::Float64 = 0.0 
    @par exps(+) for i in 1:length(X)
        d = -(X[i]-points).^2./(2*b^2)
        m = minimum(d)
        exps += m-log(b*N)+log(sum(exp(d-m)))
    end 
    return exps
end

println(kernelscore_test(128))

