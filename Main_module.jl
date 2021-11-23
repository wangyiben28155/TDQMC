module TDQMC

export Three

using Distributions
import Base.@kwdef

const L, x_num, step_t = (350.0, 1751, 13001)


@kwdef mutable struct Dynamics{T<:AbstractFloat}
    real_space::Vector{Complex{T}} = complex(collect(pdf(Truncated(Normal(0, 30), -Inf64, Inf64), LinRange(-L, L, x_num))))
    Time::Union{T,Complex{T}} = 0.0
end


@kwdef struct Parameter{T<:AbstractFloat}
    

end
 
include("Numerical_Diff_DiscreteFunc.jl")      #用来做数值微分的函数, 对离散函数(不知道解析表达式的波函数进行求解)
include("Potential.jl")
# include("")
# include("")
# include("")
# include("")



end