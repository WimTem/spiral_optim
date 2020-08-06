using LinearAlgebra, Plots, Distributions, Random

f(x, y) = begin
    -(3x + y ^ 2) * abs(sin(x) + cos(y))
end

function generate_val(dom, h)
    x = dom[1]:h:dom[2]
    y = dom[3]:h:dom[4]
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    Z = map(f, X, Y)
    return x, y, X, Y, Z
end

#Visual
dom_init = [5, 10, 4, 8]
a = generate_val(dom_init, 0.5)
contour(a[1], a[2], a[5])

function generate_rand(X, Y, m)
    rng = MersenneTwister(1234)
    randx = rand!(rng, zeros(m), X)
    randy = rand!(rng, zeros(m), Y)
    return hcat(randx, randy)'
end

function rotate(θ, x)
    return [cos(θ) -sin(θ); sin(θ) cos(θ)]*x
end

#Dom: domain
#h : stepsize
#m : amount of rand points
#θ : rotation angle
#r : radius
function SOA(dom, h, m, θ, r, n_iter)
    rng = MersenneTwister(1234)
    X, Y, Z = generate_val(dom, h)[3:5]
    X_rand = rand!(rng, zeros(m), X)
    Y_rand = rand!(rng, zeros(m), Y)
    result = hcat(X_rand, Y_rand)'#Row 1: Xval, Row2: Yval
    Xc = result[:,find_Xc(result)]
    for i = 1:n_iter
        for j = m-1:-1:0
            x = r*rotate(θ,result[:,end-j]) - r*rotate(θ, Xc) + I(2)*Xc 
            result = hcat(result,x)
        end 
        Xc = result[:, end - 5 + find_Xc(result[:, end-4:end])]
    end
    return result
end

result = SOA(dom_init, 0.1, 5, π/9, 0.80, 30)
result_x = result[1,:]
result_y = result[2,:]
contour(a[1], a[2], a[5])
scatter!(result_x[1:5:end], result_y[1:5:end])
savefig("SOA.pdf")
scatter!(result_x[2:5:end], result_y[2:5:end])
scatter!(result_x[3:5:end], result_y[3:5:end])
scatter!(result_x[4:5:end], result_y[4:5:end])
scatter!(result_x[5:5:end], result_y[5:5:end])

function find_Xc(x)
    result = 1e3
    index = []
    for i = 1:1:length(x[1,:])
        if f(x[1, i], x[2,i]) < result
            result = f(x[1, i], x[2,i])
            index = push!(index, i)
        end
    end
    return index[end]
end
