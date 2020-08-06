using LinearAlgebra, Plots, Distributions, Random, ColorSchemes

f(x, y) = begin
    (3x + y ^ 2) * abs(sin(x) + cos(y))
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

#Generate m random points in domain
function generate_rand(X, Y, m)
    rng = MersenneTwister(1234)
    randx = rand!(rng, zeros(m), X)
    randy = rand!(rng, zeros(m), Y)
    return hcat(randx, randy)'
end

#Rotate points around point with max fitness value (Xc)
function rotate(θ, x)
    return [cos(θ) -sin(θ); sin(θ) cos(θ)]*x
end

function find_Xc(x)
    result = -1e3
    index = []
    for i = 1:length(x[1,:])
        if f(x[1, i], x[2,i]) > result
            result = f(x[1, i], x[2,i])
            index = x[:,i]
        end
    end
    return index
end

#Dom: domain
#h : stepsize
#m : amount of rand points
#θ : rotation angle
#r : radius
function SOA(dom, h, m, θ, r, n_iter)
    #Generate m random points from domain, store in matrix A
    rng = MersenneTwister(1234)
    X, Y, Z = generate_val(dom, h)[3:5]
    X_rand = rand!(rng, zeros(m), X)
    Y_rand = rand!(rng, zeros(m), Y)
    A = hcat(X_rand, Y_rand)' # Xval = A[1,:], Yval = A[2,:]
    #Find Xc, critical point with highest fitness value
    Xc = find_Xc(A)

    for i = 1:n_iter
        B = Array{Float64}(undef, 2, 0)

        #Generate next m points
        B = r*rotate(θ, A[:, end - m+1 : end]) .- r*rotate(θ, Xc) .+ Xc

        #Add new points to A
        A = hcat(A, B)

        #From last m points, find new Xc
        Xc = find_Xc(B)

    end
    contour(a[1], a[2], a[5])
    scatter!(A[1,1:m:end], A[2,1:m:end])
    scatter!(A[1,2:m:end], A[2,2:m:end])
    scatter!(A[1,3:m:end], A[2,3:m:end])
    savefig("SOA.pdf")
    
    return println("Optimimum at: ", A[:,end])
end

result = SOA(dom_init, 0.1, 3, π/6, 0.87, 25)