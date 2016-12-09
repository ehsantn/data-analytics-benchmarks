
function kernelscore_test(n)
    X = rand(n)
    points = [-1.0, 2.0, 5.0]
    N = size(points,1)
    b = 0.5
    exps::Float64 = 0.0
    tic()
    for i in 1:length(X)
        d = -(X[i]-points).^2./(2*b^2)
        m = minimum(d)
        exps += m-log(b*N)+log(sum(exp(d-m)))
    end
    toc()
    return exps
end

n = 128
if length(ARGS)==1
    n = parse(Int, ARGS[1])
end

println(kernelscore_test(n))

