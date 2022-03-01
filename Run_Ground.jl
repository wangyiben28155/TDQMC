module Run_Ground

include("Main_module.jl")
using .TDQMC

P = Parameter{Float64,Int64}()
Dy = Dynamics{Float64}()

parallel_CTE!(P, Dy)

Sum_Energy = sum(Dy.Energy)/length(Dy.Energy)

println("The Total Energy is $(Sum_Energy)")
end