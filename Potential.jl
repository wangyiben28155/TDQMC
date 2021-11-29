module Potential_Matrix                                                  #这里的函数使用的

export V_Matrix, T_Matrix, V_0_Matrix

using ..TDQMC
using SparseArrays

const ω_0 = 0.148

# 导波函数部分
V_ne(x::T; α::T = 1.0) where T<:AbstractFloat = -2.0/sqrt(α+x^2)       #定义势能项, 因为这里不需要对矩阵进行运算, 就不需要用
V_ee(x_t::T; x::Vector{T},β::T = 1.0) where T<:AbstractFloat = 1.0 ./sqrt.(β .+ (x .- x_t) .^ 2 )


Envelope(t::T; T_0::T = 2pi/ω_0, ξ_0::T = 0.1) where T <:AbstractFloat = t<=3*T_0 ? ξ_0 * sin(pi*t/(6*T_0))^2 : ξ_0         
E(t::T; ϵ::Function =  Envelope , ω::T = ω_0) where T<:AbstractFloat = ϵ(t) * sin(ω*t)                   #定义电场

 
Operator(x::T,t::T) where T<:AbstractFloat = V_0(x) - x*E(t)              #自由演化的话可以不加电场


function V_Matrix(P::Parameter, Wave::Dynamics; V::Function = Operator)
    local V_operator = @. exp(-1im * P.Δt * V(P.sampling, Wave.Time))                                        #这里势能项是与时间有关的,所以需要引入Wave.Time

    return sparse(1:P.N, 1:P.N, V_operator)
end

function T_Matrix(P::Parameter; T::Function = Momentum_T)
    local T_operator = T(P)

    return sparse(1:P.N, 1:P.N, T_operator)
end

end