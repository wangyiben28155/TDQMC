module TDQMC

export Dynamics, Parameter, extend_num, stitch_Matrix, parallel_Evolution!, parallel_CTE!, 
        Group_Energy, HHG, choose

using Distributions, Random, SparseArrays                                   #用来初始化初始的波函数和系综粒子分布
import Base.@kwdef

const Ensemble_num, Electron_num = (500, 1)                                 #设定一些计算的参数
const L, x_num, step_t, Δt = (30.0, 3001, 300, 0.05 - 0.05im)
const μ, σ = (zeros(Float64, Electron_num), 1.5 * ones(Float64, Electron_num))       #这里考虑到自旋和导波函数的初始化,我们对每组系综中代表第n个电子的粒子进行轨迹的初始化的时候,导波函数应该不是相同的,否则斯莱特行列式会变成零
const spin = [1]


function initializer_guideWave(n::T; μ::AbstractVector = μ, σ::AbstractVector = σ) where {T<:Integer}                                           #这里因为算的是一维的波函数, 所以返回的矩阵为三维矩阵,第一维为波函数,第二维为系综粒子数,第三维为电子数
    return complex(sqrt.(pdf(Normal(μ[n], σ[n]), LinRange(-L, L, x_num))))
end


function generate_distribution(n::Int; μ::AbstractVector = μ, σ::AbstractVector = σ)
    Random.seed!(n)                                                       #确定下一次使用rand函数的撒点是固定的, 每次调用程序也能保证下面的分布是确定的
    Initial_Distriution = rand(Normal(μ[n], σ[n]), Ensemble_num)                #将确定的撒点分布储存起来, 用于下面实例化的赋值

    return Initial_Distriution'
end


@kwdef mutable struct Dynamics{T<:AbstractFloat}
    Trajectory::Matrix{T} = vcat(generate_distribution.(1:Electron_num)...)             #先给定一个初始的轨迹分布, 目前只是计算氢原子的情况,
    #并且之后代不同电子之间粒子的位置之间的运算会有矢量化运算, 所以将其放在列坐标加快运算 
    Guide_Wave::Matrix{<:Vector{<:Complex{T}}} = [initializer_guideWave(i) for i in 1:Electron_num, j in 1:Ensemble_num]     #这里得到的分布是概率密度的开平方作为初始的波函数
    #进行演化
    Energy::Vector{T} = zeros(T, Ensemble_num)
    Time::Vector{Union{T,Complex{T}}} = zeros(typeof(Δt), Ensemble_num)
    Displace::Array{T,3} = zeros(T, (step_t+1, Ensemble_num, Electron_num))
    Index::Vector{Vector{<:Integer}} = [Int64[0] for i in 1:Ensemble_num]                #在-P.scope到P.scope之内电子轨迹的索引
    In_num::Vector{<:Integer} = ones(Int64, Ensemble_num)                              #边界内电子的数目
end


@kwdef struct Parameter{T1<:AbstractFloat,T2<:Integer}                        #用来控制计算参数的
    electron::T2 = Electron_num
    Group::T2 = Ensemble_num
    Spin::Vector{T2} = spin
    space_N::T2 = x_num                                                       #划分的格点的总数,后面做离散傅里叶变换的时候会用得到
    scope::T1 = L                                                             #确定波函数的计算范围为-scope到+scope
    Δx::T1 = 2 * scope / (space_N - 1)                                              #波函数的离散的空间间隔
    Twice_Δx²::T1 = 2.0 * (Δx)^2
    sampling::LinRange{T1} = LinRange(-scope, scope, space_N)
    Δt::Union{T1,Complex{T1}} = Δt                                         #划分的时间间隔, 尝试时间迭代区间, 因为考虑到虚时演化, 所以类型设定为复数
    step_t::T2 = step_t
end

include("Numerical_Diff_DiscreteFunc.jl")                        #用来做数值微分的函数, 对离散函数(不知道解析表达式的波函数进行求解)
include("Stitching.jl")
include("Find_nearest.jl")
include("Potential.jl")
include("Physical_quantity.jl")
include("Crank_Nicolson.jl")
include("rejection.jl")
include("Trajectory.jl")
include("Evolution.jl")
include("Evolution_complex.jl")
include("Parallel_Calculation.jl")


using .Stitching
using .Evolution           #把包中的变量空间导入到主模块当中export出去使用
using .Trajectory
using .Quantity
using .Rejection
using .Parallelize



end