module Trajectory

export Movement

using ..TDQMC
using ..TDQMC.Discrete_Func_diff               #使用此文件中定义的数值微分的函数
using ..TDQMC.Find_nearest_k
using NumericalAnalysis: Polynomial.Lagrange   #然后对不在格点上的轨迹进行插值求得其导数等等

function find_lattice(x::Vector{T}, P::Parameter, Dy::Dynamics; k::Integer = 5) where {T<:AbstractFloat}     #这个对应的是粒子的轨迹矢量,用来得到相应的波函数部分
    local indexVec_k::Vector{<:Vector{<:Integer}} = find_k_index.(x, x = P.sampling, k = k)
    local Wave_Vector::Vector{<:Vector{<:Complex{T}}} = fill(zeros(Complex{T}, k), P.electron)

    for i = 1:P.e
        @inbounds Wave_Vector[i] = Dy.Guide_Wave[i][indexVec_k]
    end
end


function Derivative_Trajectory()
    local a

end


function Slater_determinant(P::Parameter, Dy::Dynamics)           #通过此函数得到交叉关联的波函数, 按理来说有多少个电子就应该有多少个坐标,对应Electron_num维度的电子波函数
    local symmetric_Wave



    return symmertric_Wave
end


function Movement(P::Parameter, Dy::Dynamics)                     # 这里我们使用欧拉法即可
    local 







end