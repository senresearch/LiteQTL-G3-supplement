using Statistics
function benchmark(nrep::Int64,f::Function,x...)
    ns2sec = 1.0e-9

    res = Array{Float64}(undef, nrep)

    for i=1:nrep
        start = time_ns()
        f(x...)
        finish = time_ns()
        res[i] = (finish - start)*ns2sec
    end

    return  [minimum(res) quantile(res,[0.25  0.5 0.75]) maximum(res)]

end