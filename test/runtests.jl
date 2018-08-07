using ValueFunctionSolver
using Base.Test

# write your own tests here
# test indexing 

I1 = 2
I2 = 3 
I3 = 4
X = [i1*i2*i3 for i1=1:I1, i2=1:I2, i3=1:I3]
Xlong = X[:]
for i1=1:I1, i2=1:I2, i3=1:I3
    @test X[i1, i2, i3] == Xlong[index_vcattensor([i1, i2, i3], [I1, I2, I3])]
end


xs = [linspace(0., d*.2, d+1) for d = 1:5]
ws = [ones(length(xs[d]))/length(xs[d]) for d = 1:5]
xgrid, w = tensor_grid(xs,ws)

