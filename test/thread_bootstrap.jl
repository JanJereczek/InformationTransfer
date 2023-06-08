using .Threads

function test_computation()
    X = rand(100, 100000)
    Y = zeros(size(X)...)
    for i in axes(X, 1)
        Y[i, :] .= X[i, :] .^ 2
    end
    return X, Y
end

function threaded_test_computation()
    X = rand(100, 100000)
    Y = zeros(size(X)...)
    @threads for i in axes(X, 1)
        Y[i, :] .= X[i, :] .^ 2
    end
    return X, Y
end

@btime test_computation()
@btime threaded_test_computation()

X, Y = threaded_test_computation()
Y ≈ X .^ 2