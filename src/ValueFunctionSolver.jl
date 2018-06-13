module ValueFunctionSolver
    export gss, it_solve, index_vcattensor, tensor_grid

    const invphi = (sqrt(5) - 1) / 2
    const invphi2 = (3 - sqrt(5)) / 2
    """
    Golden section search.

    Given a function f with a single local minumum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the minimum with d-c <= tol.

    example:
    f(x) = (x-2)^2
    a = 1
    b = 5
    tol = 1e-5
    bracket = gss(f, a, b, tol = tol)
    result:
    (1.9999959837979107, 2.0000050911830893)
    """
    function gss(f,a::Float64,b::Float64; tol::Float64=1e-5)
        (a,b)=(min(a,b),max(a,b))
        h = b - a
        if h <= tol
            return (a,b)
        end
        # required steps to achieve tolerance
        n = Int(ceil(log(tol/h)/log(invphi)))
        c = a + invphi2 * h
        d = a + invphi * h
        yc = f(c)
        yd = f(d)

        for k in 1:n
            if yc < yd
                b = d
                d = c
                yd = yc
                h = invphi*h
                c = a + invphi2 * h
                yc = f(c)
            else
                a = c
                c = d
                yc = yd
                h = invphi*h
                d = a + invphi * h
                yd = f(d)
            end
        end
        if yc < yd
            return (a,d)
        else
            return (c,b)
        end
    end

    function it_solve(Vguess, update_V!; tol = 1e-5, maxit = 1000, verbose = true)
        V = copy(Vguess)
        dist = tol*1e2
        for it = 1:maxit
            if verbose
                println("Iteration $it")
            end
            update_V!(V)
            dist = maximum(abs, V-Vguess)
            if verbose
                println("Distance $dist")
            end
            if isnan(dist)
                throw("Distance NaN")
            end
            if dist < tol
                return V, dist
            end
            Vguess[:] = V[:]
        end
        println("Did not reach specified tolerance within maximum number of iterations. Current dist: $dist")
        return V, dist
    end

    """
    index_vcattensor(subtensor::Array{Int}, Ntensor::Array{Int})
    
    Index into vcat of tensor grid. Subtensor position in each dimension. First dimension moves first.

    """
    function index_vcattensor(subtensor::Array{Int}, Ntensor::Array{Int})
       myind = 0
       # iterate up to second to last index
       for i in length(subtensor):-1:2
           myind += (subtensor[i] - 1)*prod(Ntensor[1:i-1])
       end
       myind += subtensor[1]
       return myind
   end

   """
   tensor_grid(xs)

   Vcat of tensor grid. First dimension moves first.
   where xs is a collection of unidimensional grids
   """
   function tensor_grid(xs)
       N = [length(xs[i]) for i = 1:length(xs)]
       totN = prod(N)
       d = length(N)
       xmat = Array{eltype(xs[1])}(totN,d)
       for i in 1:d
           repout = Int(totN/prod(N[1:i]))
           repin = Int(totN/N[i]/repout)
           xmat[:, i] = repeat(xs[i], inner=repin, outer=repout)
       end
       return xmat
end
end # module
