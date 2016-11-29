using HPAT

@acc hpat function kernelscore_test(n)
    X = rand(n)
    points = [-1.0, 2.0, 5.0]
    N = size(points,1)
    b = 0.5
    exps::Float64 = 0.0 
    @par exps(+) for i in 1:length(X)
        d = -(X[i]-points).^2./(2*b^2)
        exps += minimum(d)-log(b*N)+log(sum(exp(d-minimum(d))))
    end 
    return exps
end

println(kernelscore_test(128))

