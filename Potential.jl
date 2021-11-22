module Potential_Matrix                                                  #根据引用流和功能类比,这里相当于是厨房,用于输出最基本的用来计算需要用到的函数

export V_Matrix, T_Matrix, V_0_Matrix

using ..TDQMC
using SparseArrays

const ω_0 = 0.148

V_0(x::T; α::T = 1.0) where T<:AbstractFloat = -1.0/sqrt(α+x^2)                                             #软核势,之后势能项还需要

Envelope(t::T; T_0::T = 2pi/ω_0, ξ_0::T = 0.1) where T <:AbstractFloat = t<=3*T_0 ? ξ_0 * sin(pi*t/(6*T_0))^2 : ξ_0         
E(t::T; ϵ::Function =  Envelope , ω::T = ω_0) where T<:AbstractFloat = ϵ(t) * sin(ω*t)                      #定义电场

Operator(x::T,t::T) where T<:AbstractFloat = V_0(x) - x*E(t)


Momentum_T(P::Parameter) =@. exp(-1im * P.Δt * pi^2 * P.frequency_space^2)                                  #这里返回的是一个向量


function V_0_Matrix(P::Parameter ; V::Function = V_0)                                                       #这里其实可以把它和下面的函数合并,设定时间为零即可,但因为这样做简单直观并且没有冗余的操作,性能上应该会好一些,而且就几行占用内存也不多.
    local V_operator =@. exp(-1im * P.Δt * V(P.sampling)) 

    return sparse(1:P.N, 1:P.N, V_operator)
end

function V_Matrix(P::Parameter, Wave::wave_function; V::Function = Operator)
    local V_operator =@. exp(-1im * P.Δt * V(P.sampling, Wave.Time))                                        #这里势能项是与时间有关的,所以需要引入Wave.Time
    
    return sparse(1:P.N, 1:P.N, V_operator)
end

function T_Matrix(P::Parameter; T::Function = Momentum_T)
    local T_operator = T(P)

    return sparse(1:P.N,1:P.N,T_operator)
end

end