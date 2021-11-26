module Potential_Matrix                                                  #根据引用流和功能类比,这里相当于是厨房,用于输出最基本的用来计算需要用到的函数

export V_Matrix, T_Matrix, V_0_Matrix

using ..TDQMC
using SparseArrays

const ω_0 = 0.148

V_ne(x::T; α::T = 1.0) where T<:AbstractFloat = -2.0/sqrt(α+x^2)       #定义势能项
v_ee(x::T, ; β::T = 1.0) where T<:AbstractFloat = 1.0/sqrt(β+(x-y)^2 )


Envelope(t::T; T_0::T = 2pi/ω_0, ξ_0::T = 0.1) where T <:AbstractFloat = t<=3*T_0 ? ξ_0 * sin(pi*t/(6*T_0))^2 : ξ_0         
E(t::T; ϵ::Function =  Envelope , ω::T = ω_0) where T<:AbstractFloat = ϵ(t) * sin(ω*t)                      #定义电场



Operator(x::T,t::T) where T<:AbstractFloat = V_0(x) - x*E(t)


end