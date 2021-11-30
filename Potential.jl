module Potential_Matrix                                                  #这里的函数使用的

export V_operator

using ..TDQMC
using SparseArrays, LinearAlgebra

const ω_0 = 0.148

# 导波函数部分
V_ne(x::T; α::T = 1.0) where {T<:AbstractFloat} = -2.0 / sqrt(α + x^2)                         #定义势能项, 因为这里不需要对矩阵进行运算, 就不需要用
V_ee(x_t::T; x::Vector{T}, β::T = 1.0) where {T<:AbstractFloat} = @. 1.0 / sqrt(β + (x - x_t)^2)

Envelope(t::T; T_0::T = 2pi / ω_0, ξ_0::T = 0.1) where {T<:AbstractFloat} = t <= 3 * T_0 ? ξ_0 * sin(pi * t / (6 * T_0))^2 : ξ_0
Electric_Field(t::T; ϵ::Function = Envelope, ω::T = ω_0) where {T<:AbstractFloat} = ϵ(t) * sin(ω * t)   #定义电场


function Vne_Operator(P::Parameter) where {T<:AbstractFloat} #注意这个函数可以作为并行的基本单元,输入需要加入控制参数方便并行,注意这里输入的矢量是随时间变化的
    local Matrix_Vne::Matrix{T} = zeros(T, (P.space_N, P.Electron))

    Matrix_Vne = V_ne.(P.sampling)         #这里得到的是核势能项加上电场项对应的势能矩阵
    Matrix_Vne *= ones(T, (1, P.Electron))

    return Matrix_Vne
end


function Vee_Operator(x_t::Vector{T}, P::Parameter) where {T<:AbstractFloat}             #此函数返回的是Vector类型, 对应空间上势能的分布, x_t为不同系综粒子的轨迹,这里可以结合不同的组进行并行
    local Matrix_Vee::Matrix{T} = zeros(T, (P.space_N, P.Electron))                      #这个矩阵是用来存放势能函数的空间分布的
    local Sum_Matrix::Matrix{T} = ones(T, (P.Electron, P.Electron)) .- Diagonal(ones(T, P.Electron))

    Matrix_Vee = hcat(V_ee.(x_t, x = P.sampling)...)
    Matrix_Vee *= Sum_Matrix

    return Matrix_Vee
end

function Vtd_Operator(x_t::Vector{T}, P::Parameter, t::T; Field::Function = Electric_Field) where {T<:AbstractFloat} #随时间变化的势能项, 这里之所以拆开是为了方便后面的矩阵构造时候能更节省时间
    local 

end

function Vsum_Operator(x_t::Vector{T}, P::Parameter, t::T; Field::Function = Electric_Field) where {T<:AbstractFloat}



end


end