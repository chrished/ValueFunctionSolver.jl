# ValueFunctionSolver
Some tools typically useful for numerical solution of systems of equations

1. index_vcattensor(idim::Array{Int,1}, Idim::Array{Int,1}) 
Function to index into the vcat of a D-dimensional Array (useful for multidimensional grids).
Usage:
    ```julia
    X = [i1*i2*i3 for i1=1:I1, i2=1:I2, i3=1:I3]
    Xlong = X[:]
    X[i1, i2, i3] == Xlong[index_vcattensor([i1, i2, i3], [I1, I2, I3])]
    ```

