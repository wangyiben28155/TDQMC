module Run_Ground

include("Main_module.jl")
using .TDQMC

using Dates

tic = now()

P = Parameter{Float64,Int64}()
Dy = Dynamics{Float64,Int64}()

parallel_CTE!(P, Dy)

Sum_Energy = sum(Dy.Energy) / length(Dy.Energy)

toc = now()

println("The Total Energy is $(Sum_Energy)")

println("$(convert(DateTime, toc-tic))")

end