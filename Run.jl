include("Main_module.jl")
using .TDQMC

P = Parameter{Float64,Int64}()
Dy = Dynamics{Float64}();

parallel_CTR!(P, Dy)