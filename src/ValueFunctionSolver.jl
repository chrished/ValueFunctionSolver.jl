module ValueFunctionSolver
    export gss, it_solve

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

    function it_solve(Vguess, update_V!; tol = 1e-5, maxit = 1000)
        V = copy(Vguess)
        dist = tol*1e2
        for it = 1:maxit
            update_V!(V)
            dist = maximum(abs, V-Vguess)
            if dist < tol
                return V
            end
            Vguess[:] = V[:]
        end
        println("Did not reach specified tolerance within maximum number of iterations. Current dist: $dist")
        return V, dist
    end
end # module
