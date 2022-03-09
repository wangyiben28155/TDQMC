module Potential_Matrix                                                  #这里的函数使用的

export V_ne, V_ee, Vtd_Operator

using ..TDQMC
import ..TDQMC.Electron_num
using SparseArrays, LinearAlgebra

#import PhysicalConstants.CODATA2018: a_0, ħ, m_e, c_0

#const ω_unit = ħ / (m_e * a_0^2)
#const λ = 800E-9    # 800nm波长
#const ω_0 = round(((2pi * c_0 / λ) / ω_unit).val, digits = 3) #取小数点后3位
const λ = 800   # 800nm波长
const ω_0 = 0.148#round(45.56/λ, digits = 3)
const T_0 = 2pi / ω_0
const ξ_0 = 0.1

# 导波函数部分
V_ne(x::T; α::T = 1.0) where {T<:AbstractFloat} = -float(Electron_num) / sqrt(α + x^2)                         #定义势能项, 因为这里不需要对矩阵进行运算, 就不需要用
Envelope(t::T; ξ_0::T = ξ_0, ω::T = ω_0, T_0::T = T_0) where {T<:AbstractFloat} = t<=3T_0 ? ξ_0 * sin(ω * t / 12)^2 : ξ_0
Electric_Field(t::T; ϵ::Function = Envelope, ω::T = ω_0) where {T<:AbstractFloat} = ϵ(t) * sin(ω * t)   #定义电场

V_ee(x_t::T; x::AbstractVector{T}, β::T = 0.2) where {T<:AbstractFloat} = @. 1.0 / sqrt(β + (x - x_t)^2)

function Vee_Operator(x_t::Vector{T}, P::Parameter) where {T<:AbstractFloat}             #此函数返回的是Vector类型, 对应空间上势能的分布, x_t为不同系综粒子的轨迹,这里可以结合不同的组进行并行
    local Particle_num = length(x_t)
    local Matrix_Vee::Matrix{T} = zeros(T, (P.space_N, Particle_num))                      #这个矩阵是用来存放势能函数的空间分布的
    local Sum_Matrix::Matrix{T} = ones(T, (Particle_num, Particle_num)) .- Diagonal(ones(T, Particle_num))

    Matrix_Vee = hcat(V_ee.(x_t, x = P.sampling)...)
    Matrix_Vee *= Sum_Matrix             #因为考虑到每一次时间演化的时候都要重新计算这个矩阵, 所以我们还是通过最快速和简单的矩阵运算来得到电子之间的相互作用的势能, 进行类型转换虽然可能方便后面的代码直接使用点乘的操作, 方便代码的书写, 但是从计算耗时和底层来说, 点乘和循环并没有什么样本质的差别,反倒是进行类型的转换不仅会耗费内存空间进行深拷贝, 而且每次进行运算时候都进行类型转换很耗费时间.

    return Matrix_Vee
end


function Vtd_Operator(x_t::Vector{T}, P::Parameter, t::Union{T,Complex{T}};
    Field::Function = Electric_Field) where {T<:AbstractFloat} #随时间变化的势能项, 这里之所以拆开是为了方便后面的矩阵构造时候能更节省时间

    local Matrix_td::Matrix{T} = Vee_Operator(x_t, P)                       #多电子势能项,矩阵大小为(P.space_N, Particle_num)

    if imag(P.Δt) == 0
        Matrix_td .+= Field(t) .* P.sampling
    end
    #偶极近似下与电场相互作用的势能项,这里因为选择的电场为零,这样省去的计算的时间,之后电场不为零的情况下可以把判断的条件给去掉

    return Matrix_td                                            #返回的是随时间变化的势能矩阵,矩阵的每一列代表的是每一个电子受到的势能作用(除原子核之外)
end

end
