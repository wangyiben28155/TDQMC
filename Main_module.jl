module TDQMC

export wq

using Distributions
import Base.@kwdef

const L, x_num, step_t = (350.0, 1751, 13001)


@kwdef mutable struct Dynamics{T<:AbstractFloat}
    real_space::Vector{Complex{T}} = complex(collect(pdf(Truncated(Normal(0, 30), -Inf64, Inf64), LinRange(-L, L, x_num))))
    Time::Union{T,Complex{T}} = 0.0
end


@kwdef struct Parameter{T<:AbstractFloat}
    

end
 
include("Numerical_Diff_DiscreteFunc.jl")
include("Potential.jl")
include("")
include("")
include("")
include("")



end