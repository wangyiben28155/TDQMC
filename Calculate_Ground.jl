include("Main_module.jl")
using .TDQMC
using PyPlot


P = Parameter()
D = Dynamics()
Caculation!(P, D)

